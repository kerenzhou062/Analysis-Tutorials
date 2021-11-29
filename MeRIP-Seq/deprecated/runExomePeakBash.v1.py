#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict
from itertools import combinations

#usage: runExomePeakBash.py or runExomePeakBash.py <bam dir>

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-gtf', action='store', type=str, required=True,
                    help='gtf annotation gtf')
parser.add_argument('-grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('-grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input bam directory (HepG2_shWTAP_IP_rep1.sorted.bam)')
parser.add_argument('-log', action='store', type=str, required=True,
                    help='log directory')
parser.add_argument('-memory', action='store', type=str,
                    default='50G', help='memory used for sbatch')
parser.add_argument('-allRep', action='store_true',
                    default=False, help='run with only on all relicates')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output folder')
parser.add_argument('-public', action='store', type=str,
                    help='exp prefix of public input (HepG2_WT_input)')
parser.add_argument('-script', action='store', type=str, required=True,
                    help='directory for storing scripts')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

##public arguments
memory = args.memory
gtf = os.path.realpath(args.gtf)
bashDir = os.path.realpath(args.script)
logDir = os.path.realpath(args.log)
outputDir = os.path.realpath(args.output)
inputDir = os.path.realpath(args.input)

if bool(args.grepKept):
    kept = True
    regex = re.compile(r'{0}'.format(args.grepKept))
elif bool(args.grepExpel):
    kept = False
    regex = re.compile(r'{0}'.format(args.grepExpel))
else:
    regex = False

#filename: HepG2_control_input_rep1.sorted.bam, HepG2_control_input_rep2.sorted.bam
bamFiles = sorted(glob(os.path.join(inputDir, '**', '*.bam'), recursive=True))
bamList = list()
for bamFile in bamFiles:
    fileName = os.path.split(bamFile)[-1]
    if bool(regex):
        if kept:
            if bool(regex.search(fileName)) is False:
                continue
        else:
            if bool(regex.search(fileName)) is True:
                continue
    bamList.append(bamFile)

publicFlag = False
bamDict = defaultdict(dict)
publicDict = defaultdict(dict)
for bam in bamList:
    basename = os.path.basename(bam)
    cleanBasename = basename.split('.')[0]
    #HepG2 control IP rep1
    basenameList = cleanBasename.split('_')
    ip = basenameList[2]
    rep = basenameList[3]
    #HepG2_control
    baseExpName = '_'.join(basenameList[0:2])
    if bool(args.public):
        inputName = '_'.join(basenameList[0:3])
        if inputName == args.public:
            publicDict[rep] = bam
            publicFlag = True
        else:
            if baseExpName not in bamDict:
                bamDict[baseExpName] = defaultdict(dict)
            bamDict[baseExpName][ip][rep] = bam
    else:
        if baseExpName not in bamDict:
            bamDict[baseExpName] = defaultdict(dict)
        bamDict[baseExpName][ip][rep] = bam

if publicFlag is False and bool(args.public) is True:
    print("No public input found!")
    sys.exit(0)

eachRTemplate = '''#!/usr/bin/env Rscript
#exomepeak Script
#R script
#Define parameters and load library
library("exomePeak")

setwd("{outputDir}")

GTF_ANNO = "{gtf}"
# peak calling
exomepeak(GENE_ANNO_GTF=GTF_ANNO, 
  IP_BAM=c({ipBam}),
  INPUT_BAM=c({inputBam}), 
  EXPERIMENT_NAME="./{rScriptBasename}")

sessionInfo()
'''

sbatchTemplate = '''#!/bin/sh
#SBATCH --job-name={rScriptBasename}    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail  
#SBATCH -n 1                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem={memory}                      # Amount of memory in GB
#SBATCH --time=120:00:00               # Time limit hrs:min:sec
#SBATCH --output={logDir}/{rScriptBasename}.log   # Standard output and error log

BASE={bashDir}

Rscript $BASE/{rScriptBasename}.r
'''

runBashScript = os.path.join(bashDir, "runExomePeak.sbatch")
with open(runBashScript, 'w') as runBashScriptO:
    runBashScriptO.write('#!/bin/sh\n')
    runBashScriptO.write('BASE={0}\n\n'.format(bashDir))
    for baseExpName in sorted(bamDict.keys()):
        ipRepList = sorted(bamDict[baseExpName]['IP'].keys())
        if publicFlag is True:
            inputRepList = sorted(publicDict.keys())
        else:
            inputRepList = sorted(bamDict[baseExpName]['input'].keys())
        if len(ipRepList) == len(inputRepList):
            repDict = defaultdict(dict)
            for rep in ipRepList:
                repNum = int(re.findall(r'rep(\w+)', rep)[0])
                repDict[repNum]['IP'] = bamDict[baseExpName]['IP'][rep]
                if publicFlag is True:
                    repDict[repNum]['input'] = publicDict[rep]
                else:
                    repDict[repNum]['input'] = bamDict[baseExpName]['input'][rep]
            repNumList = sorted(repDict.keys())
            combinList = list()
            #generate combinations: [('1',), ('2',), ('1', '2')]
            if args.allRep:
                combinList.append(set(repNumList))
            else:
                for i in list(range(1, len(repNumList)+1)):
                    combinList.extend(list(combinations(repNumList, i)))
            for combin in combinList:
                combin = list(combin)
                combinName = '-'.join(map(str, combin))
                #HepG2_control_exomePeak_rep_1
                rScriptBasename = '_'.join([baseExpName, 'rep' + combinName])
                rScript = os.path.join(bashDir, rScriptBasename+'.r')
                sbatchScript = os.path.join(bashDir, rScriptBasename + '.sbatch')
                with open (rScript, 'w') as rScriptO, open (sbatchScript, 'w') as sbatchScriptO:
                    ipBam = ','.join(list(map(lambda x:'"' + repDict[x]['IP'] + '"', combin)))
                    inputBam = ','.join(list(map(lambda x:'"' + repDict[x]['input'] + '"', combin)))
                    rScriptCont = eachRTemplate.format(**vars())
                    rScriptO.write(rScriptCont)
                    sbatchScriptCont = sbatchTemplate.format(**vars())
                    sbatchScriptO.write(sbatchScriptCont)
                runBashScriptO.write("sbatch $BASE/{rScriptBasename}.sbatch\n".format(**vars()))
        else:
            rScriptBasename = '_'.join([baseExpName, 'exomePeak_reps'])
            rScript = os.path.join(bashDir, rScriptBasename+'.r')
            sbatchScript = os.path.join(bashDir, '_'.join([baseExpName, 'exomePeak_reps.sbatch']))
            with open (rScript, 'w') as rScriptO, open (sbatchScript, 'w') as sbatchScriptO:
                ipBam = ','.join(list(map(lambda x:'"{0}"'.format(x), bamDict[baseExpName]['IP'])))
                inputBam = ','.join(list(map(lambda x:'"{0}"'.format(x), bamDict[baseExpName]['input'])))
                rScriptCont = eachRTemplate.format(**vars())
                rScriptO.write(rScriptCont)
                sbatchScriptCont = sbatchTemplate.format(**vars())
                sbatchScriptO.write(sbatchScriptCont)
            runBashScriptO.write("sbatch $BASE/{rScriptBasename}.sbatch\n".format(**vars()))
os.chmod(runBashScript, 0o755)
