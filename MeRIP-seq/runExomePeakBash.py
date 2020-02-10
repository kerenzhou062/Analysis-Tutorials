#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict
from itertools import combinations
from PubAlbum import Anno

#usage: runExomePeakBash.py or runExomePeakBash.py <bam dir>

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-gtf', action='store', type=str,
                    default='hg38v31', help='gtf annotation build or file (eg. hg38v31)')
parser.add_argument('-input', action='store', type=str,
                    default='./', help='input bam directory (HepG2_shWTAP_IP_rep1.sorted.bam)')
parser.add_argument('-index', action='store', type=str,
                    default='hg38', help='bowtie2 genome index prefix or file (eg. hg38)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    print("Running runExomePeakBash.py with defaultdict parameters...")

##public arguments
anno = Anno()
gtf = anno.gtf(args.gtf)
genomeIndex = anno.gindex(args.index)
basepath = os.path.realpath(args.input)

mainAlignDir = os.path.join(basepath, "alignment")
exomePeakDir = os.path.join(basepath, "exomePeak")
#filename: HepG2_control_input_rep1_run.sorted.bam, HepG2_control_input_rep2_run.sorted.bam
bamList = sorted(glob(os.path.join(mainAlignDir, "*.bam")))
bashDir = os.path.join(basepath, 'runBash')
logDir = os.path.join(bashDir, 'log')

bamDict = defaultdict(dict)
for bam in bamList:
    basename = os.path.basename(bam)
    cleanBasename = re.sub(r'\..+$', '', basename)
    #HepG2 control IP rep1 run
    basenameList = cleanBasename.split('_')
    ip = basenameList[2]
    rep = "_".join(basenameList[3:5])
    #HepG2_control
    baseExpName = '_'.join(basenameList[0:2])
    if baseExpName not in bamDict:
        bamDict[baseExpName] = defaultdict(dict)
        bamDict[baseExpName][ip][rep] = basename
    else:
        bamDict[baseExpName][ip][rep] = basename

eachRTemplate = '''#!/usr/bin/env Rscript
#exomepeak Script
#R script
#Define parameters and load library
library("exomePeak")

base = "{basepath}"
setwd(paste(base, "alignment", sep = "/"))

GTF_ANNO = "{gtf}"
# peak calling
exomepeak(GENE_ANNO_GTF=GTF_ANNO, 
  IP_BAM=c({ipBam}),
  INPUT_BAM=c({inputBam}), 
  EXPERIMENT_NAME=paste("../exomePeak", "{rScriptBasename}", sep = "/"))
'''

sbatchTemplate = '''#!/bin/sh
#SBATCH --job-name={rScriptBasename}    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail  
#SBATCH -n 1                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=100G                      # Amount of memory in GB
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output={bashDir}/log/{rScriptBasename}.log   # Standard output and error log

base={bashDir}

Rscript $base/{rScriptBasename}.r
'''

runBashScript = os.path.join(bashDir, "runExomePeak.sh")
with open(runBashScript, 'w') as runBashScriptO:
    runBashScriptO.write('#!/bin/sh\n')
    runBashScriptO.write('base={0}\n\n'.format(bashDir))
    for baseExpName in sorted(bamDict.keys()):
        IPRepList = sorted(bamDict[baseExpName]['IP'].keys())
        inputRepList = sorted(bamDict[baseExpName]['input'].keys())
        if len(IPRepList) == len(inputRepList):
            repDict = defaultdict(dict)
            for rep in IPRepList:
                repNum = re.findall(r'rep(\w+)_', rep)[0]
                repDict[repNum]['IP'] = bamDict[baseExpName]['IP'][rep]
                repDict[repNum]['input'] = bamDict[baseExpName]['input'][rep]
            repNumList = sorted(repDict.keys())
            combinList = list()
            #generate combinations: [('1',), ('2',), ('1', '2')]
            for i in list(range(1, len(repNumList)+1)):
                combinList.extend(list(combinations(repNumList, i)))
            for combin in combinList:
                combin = list(combin)
                combinName = '-'.join(combin)
                #HepG2_control_exomePeak_rep_1
                rScriptBasename = '_'.join([baseExpName, 'exomePeak_rep_' + combinName])
                rScript = os.path.join(bashDir, rScriptBasename+'.r')
                sbatchScript = os.path.join(bashDir, rScriptBasename + '.sbatch')
                with open (rScript, 'w') as rScriptO, open (sbatchScript, 'w') as sbatchScriptO:
                    ipBam = ','.join(list(map(lambda x:'"' + repDict[x]['IP'] + '"', combin)))
                    inputBam = ','.join(list(map(lambda x:'"' + repDict[x]['input'] + '"', combin)))
                    rScriptCont = eachRTemplate.format(**vars())
                    rScriptO.write(rScriptCont)
                    sbatchScriptCont = sbatchTemplate.format(**vars())
                    sbatchScriptO.write(sbatchScriptCont)
                runBashScriptO.write("sbatch $base/{rScriptBasename}.sbatch\n".format(**vars()))
        else:
            rScriptBasename = '_'.join([baseExpName, 'exomePeak_reps'])
            rScript = os.path.join([bashDir, rScriptBasename+'.r'])
            sbatchScript = os.path.join(bashDir, '_'.join([baseExpName, 'exomePeak_reps.sbatch']))
            with open (rScript, 'w') as rScriptO, open (sbatchScript, 'w') as sbatchScriptO:
                ipBam = ','.join(list(map(lambda x:'"{0}"'.format(x), bamDict[baseExpName]['IP'])))
                inputBam = ','.join(list(map(lambda x:'"{0}"'.format(x), bamDict[baseExpName]['input'])))
                rScriptCont = eachRTemplate.format(**vars())
                rScriptO.write(rScriptCont)
                sbatchScriptCont = sbatchTemplate.format(**vars())
                sbatchScriptO.write(sbatchScriptCont)
            runBashScriptO.write("sbatch $base/{rScriptBasename}.sbatch\n".format(**vars()))
os.chmod(runBashScript, 0o755)
