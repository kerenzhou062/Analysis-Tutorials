#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict
import datetime
from PubAlbum import Anno

#usage: runM6aTophatAlignBash.py or runM6aTophatAlignBash.py <fastq dir>

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-cpu', action='store', type=str,
                    default='10', help='threads used for tophat')
parser.add_argument('-gtf', action='store', type=str,
                    default='hg38v31', help='gtf annotation build or file (eg. hg38v31)')
parser.add_argument('-input', action='store', type=str,
                    default='./', help='input fastq directory (eg. HepG2_shWTAP_IP_rep1_run1_2.fastq)')
parser.add_argument('-index', action='store', type=str,
                    default='hg38', help='bowtie2 genome index prefix or file (eg. hg38)')
parser.add_argument('-memory', action='store', type=str,
                    default='50G', help='memory used for sbatch')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    print("Running runM6aTophatAlignBash.py with defaultdict parameters...")

#public arguments
threadNum = args.cpu
memory = args.memory
Anno = Anno()
gtf = Anno.gtf(args.gtf)
genomeIndex = Anno.gindex(args.index)
basepath = os.path.realpath(args.input)

#filename: HepG2_shWTAP_IP_rep1_run1_1.fastq, HepG2_shWTAP_IP_rep1_run1_2.fastq, HepG2_shWTAP_IP_run_rep1_NA.fastq
extsList = ["*.fastq", '*.fq', '*.fastq.gz', '*.fq.gz']
fastqList = sorted([f for ext in extsList for f in glob(os.path.join(basepath, ext))])
bashDir = os.path.join(basepath, 'runBash')
logDir = os.path.join(bashDir, 'log')
os.makedirs(logDir, exist_ok=True)

fqApp = '' # .fq.gz|.fq|.fastq|.fastq.gz
expFqDict = defaultdict(dict)
# expFqDict: exp->IP->rep->[fq1,fq2]
for fastq in fastqList:
    basename = os.path.basename(fastq)
    fqApp = re.findall(r'\..+$', basename)[0]
    cleanBasename = re.sub(r'\..+$', '', basename)
    #HepG2 shWTAP IP rep1 1
    basenameList = cleanBasename.split('_')
    #HepG2_shWTAP
    exp = '_'.join(basenameList[0:2])
    ip = basenameList[2]
    rep = basenameList[3]
    run = basenameList[4]
    pairedNum = basenameList[5]
    if ip in expFqDict[exp]:
        if pairedNum in expFqDict[exp][ip][rep]:
            expFqDict[exp][ip][rep][pairedNum].append(basename)
        else:
            expFqDict[exp][ip][rep][pairedNum] = [basename]
    else:
        expFqDict[exp][ip] = defaultdict(dict)
        expFqDict[exp][ip][rep][pairedNum] = [basename]
mainAlignDir = os.path.join(basepath, 'alignment')

eachBashContTemplate = '''#!/bin/sh
#SBATCH --job-name=tophat_{baseExpName}    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail  
#SBATCH -n {threadNum}                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem={memory}                      # Amount of memory in GB
#SBATCH --time=72:10:00               # Time limit hrs:min:sec
#SBATCH --output={basepath}/runBash/log/{baseExpName}.log   # Standard output and error log

BASE={basepath}
ALIGN_DIR={mainAlignDir}
ANNO={gtf}
GENOME_INDEX={genomeIndex}

{catCommand}

echo "align reads with tophat..."

tophat -p {threadNum} -g 1 --library-type=fr-firststrand -G $ANNO \\
  -o $ALIGN_DIR/{baseExpName} \\
  $GENOME_INDEX \\
  {fastqParam}

echo "sort bam with samtools..."

mv $ALIGN_DIR/{baseExpName}/accepted_hits.bam $ALIGN_DIR/{baseExpName}.bam
samtools sort --threads {threadNum} -m 20G -O bam \\
  -o $ALIGN_DIR/{baseExpName}.sorted.bam $ALIGN_DIR/{baseExpName}.bam
rm $ALIGN_DIR/{baseExpName}.bam
samtools index -b $ALIGN_DIR/{baseExpName}.sorted.bam

echo "All jobs done!"
'''

runBashScript = os.path.join(bashDir, "runTophatSortRenameBam.sbatch")
with open(runBashScript, 'w') as bashO:
    bashO.write('#!/bin/sh\n')
    bashO.write('base={bashDir}\n\n'.format(**vars()))
    # expFqDict: exp->IP->rep->[fq1,fq2]
    for exp in sorted(expFqDict.keys()):
        for ip in sorted(expFqDict[exp].keys()):
            ipDict = expFqDict[exp][ip]
            for rep in sorted(ipDict.keys()):
                pairedNumList = sorted(ipDict[rep].keys())
                baseExpName = '_'.join([exp, ip, rep])
                baseAlignDir = os.path.join(mainAlignDir, baseExpName)
                os.makedirs(baseAlignDir, exist_ok=True)
                ##cat fastq files for multiple runs
                pairedNum1 = pairedNumList[0]
                runFileList = sorted(ipDict[rep][pairedNum1])
                catCommand = ''
                if len(pairedNumList) > 1:
                    if len(runFileList) > 1:
                        catCommand = 'prefix="${{alignDir}}/{baseExpName}/{baseExpName}"\n'.format(**vars())
                        catFastqR1 = '${prefix}' + '_1.{fqApp}'.format(**vars())
                        catFastqR2 = '${prefix}' + '_2.{fqApp}'.format(**vars())
                        fastqR1Runs = ' '.join(list(map(lambda x: '$BASE/'+x, sorted(ipDict[rep]['1']))))
                        fastqR2Runs = ' '.join(list(map(lambda x: '$BASE/'+x, sorted(ipDict[rep]['2']))))
                        catCommand += 'echo "cat fastqs..."\ncat {fastqR1Runs} | \\ >{catFastqR1}\n\
                            cat {fastqR2Runs} | \\ > {catFastqR2}\n'.format(**vars())
                        ipDict[rep]['1'] = catFastqR1
                        ipDict[rep]['2'] = catFastqR2
                    else:
                        ipDict[rep]['1'] = '$BASE/' + ipDict[rep]['1'][0]
                        ipDict[rep]['2'] = '$BASE/' + ipDict[rep]['2'][0]
                else:
                    if len(runFileList) > 1:
                        catCommand = 'prefix="${{alignDir}}/{baseExpName}/{baseExpName}"\n'.format(**vars())
                        catFastq =  '${prefix}' + '.fastq'
                        fastqRuns = ' '.join(list(map(lambda x: '$BASE/'+x, sorted(ipDict[rep][pairedNum1]))))
                        catCommand += 'echo "cat fastqs..."\ncat {fastqRuns} > {catFastq}\n'.format(**vars())
                        ipDict[rep][pairedNum1] = catFastq
                    else:
                        ipDict[rep][pairedNum1] = '$BASE/' + ipDict[rep][pairedNum1][0]
                eachrunBashScript = os.path.join(bashDir, baseExpName + '.sh')
                with open (eachrunBashScript, 'w') as eachBashO:
                    if len(ipDict[rep]) > 1:
                        fastqParam = ' '.join(map(lambda x:'$BASE/'+ipDict[rep][x], sorted(ipDict[rep].keys())))
                    else:
                        fastqParam = '$BASE/' + ipDict[rep][pairedNum1]
                    eachBashCont = eachBashContTemplate.format(**vars())
                    eachBashO.write(eachBashCont)
                os.chmod(eachrunBashScript, 0o755)
                logFileName = '_'.join(['tophatRun', baseExpName, '.log'])
                bashO.write('sbatch $BASE/{baseExpName}.sh\n'.format(**vars()))
os.chmod(runBashScript, 0o755)
