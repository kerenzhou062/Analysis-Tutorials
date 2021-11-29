#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict
from copy import copy
import datetime
from PubAlbum import Anno

#usage: runChipPeakCallBash.py or runChipPeakCallBash.py <fastq dir>

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-blacklist', action='store', type=str,
                    default='hg38', help='ENCODE blacklist: version or file')
parser.add_argument('-chrsize', action='store', type=str,
                    default='hg38', help='chromosome size file (eg. hg38.chrom.sizes)')
parser.add_argument('-control', action='store', type=str,
                    help='The prefix name of control input sample (like:HepG2_input) \
                    (otherwise infer it automatically)')
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='threads used for \
                    callMacs2Parallel.py|callEncodeIdrParallel.py')
parser.add_argument('-extsize', action='store', type=int,
                    help='--extsize parameter (macs2, histone:147)')
parser.add_argument('-gsize', action='store', type=str,
                    default='hs', help='-g parameter (macs2, mm|hs|dm|ce)')
parser.add_argument('-input', action='store', type=str,
                    default='./', help='input fastq directory \
                    (eg. HepG2_CTCF_control_IP_rep1_run1_1.fastq)')
parser.add_argument('-idrThresh', action='store', type=float,
                    default=0.05, help='--soft-idr-threshold parameter (idr)')
parser.add_argument('-grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('-grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('-alignment', action='store', type=str,
                    default='alignment',
                    help='alignment result directory')
parser.add_argument('-bash', action='store', type=str,
                    default='runBash',
                    help='runBash directory')
parser.add_argument('-log', action='store', type=str,
                    default='log',
                    help='log directory of runBash')
parser.add_argument('-peak', action='store', type=str,
                    default='peak',
                    help='peak results directory')
parser.add_argument('-maxPeak', action='store', type=str,
                    default='300000', help='-maxPeak parameter (SPP)')
parser.add_argument('-memory', action='store', type=str,
                    default='100G', help='memory used for sbatch')
parser.add_argument('-noIDR', action='store_true',
                    default=False, help='disable idr calling')
parser.add_argument('-part', action='store', type=str,
                    default='all', help='partition of slurm server')
parser.add_argument('-pval', action='store', type=str,
                    default='1e-2', help='-p parameter (macs2)')
parser.add_argument('-rank', action='store', type=str,
                    choices=['p.value','q.value','signal.value'],
                    default='p.value', help='--rank parameter (idr)')
parser.add_argument('-shift', action='store', type=str,
                    default='0', help='--shift parameter (macs2)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

# public arguments
threadNum = args.cpu
maxPeak = args.maxPeak
gsize = args.gsize
idrThresh = args.idrThresh
memory = args.memory
partition = args.part
pval = args.pval
rank = args.rank
shift = args.shift
anno = Anno()
blacklist = anno.blacklist(args.blacklist)
chrsize = anno.gsize(args.chrsize)
basepath = os.path.realpath(args.input)
tagAlignFinalApp = '.filt.srt.nodup.final.tagAlign.gz'
tagAlignPr1App = '.filt.nodup.pr1.tagAlign.gz'
tagAlignPr2App = '.filt.nodup.pr2.tagAlign.gz'
tagAlignPoolPr1App = '.pooled.pr1.tagAlign.gz'
tagAlignPoolPr2App = '.pooled.pr2.tagAlign.gz'
tagAlignPoolApp = '.pooled.tagAlign.gz'

# file pattern: HepG2_CTCF_knockdown_IP_rep1_run1_1.fastq, HepG2_CTCF_knockdown_IP_rep1_run1_2.fastq
##public arguments
if bool(args.grepKept):
    kept = True
    regex = re.compile(r'{0}'.format(args.grepKept))
elif bool(args.grepExpel):
    kept = False
    regex = re.compile(r'{0}'.format(args.grepExpel))
else:
    regex = False

extsList = ["*.fastq", '*.fq', '*.fastq.gz', '*.fq.gz']
fastqRawList = sorted([f for ext in extsList for f in glob(os.path.join(basepath, ext))])

fastqList = list()
if bool(regex):
    for fastq in fastqRawList:
        basename = os.path.basename(fastq)
        if kept:
            if bool(regex.search(basename)) is True:
                fastqList.append(fastq)
        else:
            if bool(regex.search(basename)) is False:
                fastqList.append(fastq)
else:
    fastqList = fastqRawList

if args.alignment == 'alignment':
    mainAlignDir = os.path.join(basepath, 'alignment')
else:
    mainAlignDir = os.path.realpath(args.alignment)

if args.peak == 'peak':
    mainPeakDir = os.path.join(basepath, 'peak')
else:
    mainPeakDir = os.path.realpath(args.peak)

if args.bash == 'runBash':
    bashDir = os.path.join(basepath, 'runBash')
else:
    bashDir = os.path.realpath(args.bash)

if args.log == 'log':
    logDir = os.path.join(bashDir, 'log')
else:
    logDir = os.path.realpath(args.log)

os.makedirs(logDir, exist_ok=True)
os.makedirs(mainPeakDir, exist_ok=True)

# exp CC_SCORES dict: store extsize for macs2 (fragment size)
# read from 3col (estFragLen) CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
expInfoDict = defaultdict(dict)
expDict = defaultdict(dict)
# expDict: exp->IP->[rep1, rep2]
for fastq in fastqList:
    basename = os.path.basename(fastq)
    cleanBasename = re.sub(r'\..+$', '', basename)
    #HepG2 CTCF knockdown IP rep1 run1 1
    basenameList = cleanBasename.split('_')
    #HepG2_CTCF_knockdown
    exp = '_'.join(basenameList[0:3])
    ip = basenameList[3]
    rep = basenameList[4]
    if ip in expDict[exp]:
        if rep not in expDict[exp][ip]:
            expDict[exp][ip].append(rep)
    else:
        expDict[exp][ip] = [rep]
    # get fragment size (estFragLen) from cross-correlation ScoreFile
    if bool(args.extsize) is False:
        if exp.find('input') < 0:
            if ip == 'IP':
                ccScoreFile = '_'.join([exp, ip, rep]) + '.trim.filt.sample.15.tagAlign.gz.cc.qc'
                ccScoreFilePath = os.path.join(mainAlignDir, exp, 'IP', rep, ccScoreFile)
                with open(ccScoreFilePath, 'r') as f:
                    line = f.readline()
                    row = line.strip().split('\t')
                    estFragLen = int(row[2])
                    if 'extsize' in expInfoDict[exp]:
                        expInfoDict[exp]['extsize'] += estFragLen
                        expInfoDict[exp]['repNum'] += 1
                    else:
                        expInfoDict[exp]['extsize'] = estFragLen
                        expInfoDict[exp]['repNum'] = 1

# infering public control input samples
pubConDict = defaultdict(list)
expList = sorted(expDict.keys())
if args.control:
    #like HepG2_input
    exp = args.control
    pubConDict['basename'] = exp
    pubConNum = len(expDict[exp]['input'])
    for rep in sorted(expDict[exp]['input']):
        repPrefix = '_'.join([exp, 'input', rep])
        pubConDict['tagAlignFinal'].append(repPrefix + tagAlignFinalApp)
        pubConDict['repPr1'].append(repPrefix + tagAlignPr1App)
        pubConDict['repPr2'].append(repPrefix + tagAlignPr2App)
    pubConDict['pooled'].append(exp + '_input' + tagAlignPoolApp)
    pubConDict['pooled'].append(exp + '_input' + tagAlignPoolPr1App)
    pubConDict['pooled'].append(exp + '_input' + tagAlignPoolPr2App)
else:
    for exp in expList:
        #like HepG2_input
        if exp.find('input') > 0:
            pubConDict['basename'] = exp
            pubConNum = len(expDict[exp]['input'])
            for rep in sorted(expDict[exp]['input']):
                repPrefix = '_'.join([exp, 'input', rep])
                pubConDict['tagAlignFinal'].append(repPrefix + tagAlignFinalApp)
                pubConDict['repPr1'].append(repPrefix + tagAlignPr1App)
                pubConDict['repPr2'].append(repPrefix + tagAlignPr2App)
            pubConDict['pooled'].append(exp + '_input' + tagAlignPoolApp)
            pubConDict['pooled'].append(exp + '_input' + tagAlignPoolPr1App)
            pubConDict['pooled'].append(exp + '_input' + tagAlignPoolPr2App)
            break

expLinkDict = defaultdict(dict)

for exp in expList:
    if exp.find('input') > 0 or 'IP' not in expDict[exp]:
        continue
    if exp not in expLinkDict:
        expLinkDict[exp] = defaultdict(list)
    ipList = sorted(expDict[exp].keys())
    if len(ipList) == 1 and ipList[0] == 'input':
        continue
    # get information from IP
    repList = sorted(expDict[exp]['IP'])
    repNum = len(repList)
    expLinkDict[exp]['name'].extend(list(map(lambda x:'_'.join([exp, x, 'final']), repList)))
    expLinkDict[exp]['name'].extend(list(map(lambda x:'_'.join([exp, x, 'pr1']), repList)))
    expLinkDict[exp]['name'].extend(list(map(lambda x:'_'.join([exp, x, 'pr2']), repList)))
    expLinkDict[exp]['name'].append('_'.join([exp, 'pooled']))
    expLinkDict[exp]['name'].append('_'.join([exp, 'pr1', 'pooled']))
    expLinkDict[exp]['name'].append('_'.join([exp, 'pr2', 'pooled']))
    expLinkDict[exp]['other'] = ["--keep-dup all -B --SPMR" for i in range(len(expLinkDict[exp]['name'])) ]
    outputList = list(map(lambda x:os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', \
        '_'.join([exp, x, 'final'])), repList))
    outputList.extend(list(map(lambda x:os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', \
        '_'.join([exp, x, 'pr1'])), repList)))
    outputList.extend(list(map(lambda x:os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', \
        '_'.join([exp, x, 'pr2'])), repList)))
    outputList.append(os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', exp+'_pooled'))
    outputList.append(os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', exp+'_pr1_pooled'))
    outputList.append(os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', exp+'_pr2_pooled'))
    expLinkDict[exp]['output'] = outputList
    conExp = '_'.join([exp, 'input'])
    for ip in ipList:
        if ip == 'input':
            inputRepList = expDict[exp]['input']
            if len(inputRepList) < repNum:
                extraNum = repNum - len(inputRepList)
                extraList = [inputRepList[i] for i in range(extraNum)]
                inputRepList.extend(extraList)
                repList = inputRepList
        repList = sorted(expDict[exp][ip])
        # order: final, pr1, pr2, pooled
        ipTagAlignList = list(map(lambda x:os.path.join( "${MAIN_ALIGN_DIR}", exp, ip, \
            '_'.join([exp, ip, x]) + tagAlignFinalApp), repList))
        ipTagAlignList.extend(list(map(lambda x:os.path.join( "${MAIN_ALIGN_DIR}", exp, ip, \
            '_'.join([exp, ip, x]) + tagAlignPr1App), repList)))
        ipTagAlignList.extend(list(map(lambda x:os.path.join( "${MAIN_ALIGN_DIR}", exp, ip, \
            '_'.join([exp, ip, x]) + tagAlignPr2App), repList)))
        ipTagAlignList.append(os.path.join( "${MAIN_ALIGN_DIR}", exp, ip, exp+'_'+ip+tagAlignPoolApp))
        ipTagAlignList.append(os.path.join( "${MAIN_ALIGN_DIR}", exp, ip, exp+'_'+ip+tagAlignPoolPr1App))
        ipTagAlignList.append(os.path.join( "${MAIN_ALIGN_DIR}", exp, ip, exp+'_'+ip+tagAlignPoolPr2App))
        expLinkDict[exp][ip] = ipTagAlignList
    if len(ipList) == 1:
        conExp = pubConDict['basename']
        if pubConNum < repNum:
            extraNum = repNum - pubConNum
            tagAlignFinalList = copy(pubConDict['tagAlignFinal'])
            tagAlignFinalList.extend([tagAlignFinalList[i] for i in range(extraNum)])
            repPr1List = copy(pubConDict['repPr1'])
            repPr1List.extend([repPr1List[i] for i in range(extraNum)])
            repPr2List = copy(pubConDict['repPr2'])
            repPr2List.extend([repPr2List[i] for i in range(extraNum)])
        else:
            tagAlignFinalList = copy(pubConDict['tagAlignFinal'])
            repPr1List = copy(pubConDict['repPr1'])
            repPr2List = copy(pubConDict['repPr2'])
        inputList = list(map(lambda x:os.path.join( "${MAIN_ALIGN_DIR}", 
            conExp, 'input', x ), tagAlignFinalList))
        inputList.extend(list(map(lambda x:os.path.join( "${MAIN_ALIGN_DIR}", 
            conExp, 'input', x ), repPr1List)))
        inputList.extend(list(map(lambda x:os.path.join( "${MAIN_ALIGN_DIR}", 
            conExp, 'input', x ), repPr2List)))
        inputList.extend(list(map(lambda x:os.path.join( "${MAIN_ALIGN_DIR}", 
            conExp, 'input', x ), pubConDict['pooled'])))
        expLinkDict[exp]['input'] = inputList

# ================================
#sbatch
# ================================

sbatchMacs2IdrTemplate = '''#!/bin/bash
#SBATCH --job-name=Macs2_PeakCall_{baseExpName}    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Mail user
#SBATCH -n {threadNum}                # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p {partition}                # default queue is all if you don't specify
#SBATCH --mem={memory}                # Amount of memory in GB
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output={logPath}            # Standard output and error log
### this pipeline mainly based on ENCODE3 pipeline:
###https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c
# ================================
# load required module
# ================================
echo -e "load modules..."
module load bedtools/2.29.0
module load Python/3.6.5
module load Anaconda2/4.2.0
# ===========step 0a==============
# public bariables
# ================================
echo -e "\\nPreparing basic variables..."
BASE="{basepath}"
MAIN_ALIGN_DIR="{mainAlignDir}"
MAIN_PEAK_DIR="{mainPeakDir}"
#separated by " "
THREADS={threadNum}
GENOMESIZE="{gsize}"
IP_CHIP="{ipChipStr}"
CONTROL_CHIP="{conChipStr}"
NAME="{nameStr}"
OTHER="{otherStr}"
OUTPUT="{outputStr}"
EXTSIZE="{extsize}"
SHIFT="{shift}"
CHRSIZE="{chrsize}"
RANK="{rank}"
echo "Change wd directory to ${{MAIN_PEAK_DIR}}"
cd ${{MAIN_PEAK_DIR}}
# =========call peak with macs2 ===========
# Parallel running
# =======================
echo "Peak-calling starting & sorting peak with {rank}..."
callMacs2Parallel.py -cpu ${{THREADS}} \\
  -chrsize ${{CHRSIZE}}\\
  -ip ${{IP_CHIP}} \\
  -input ${{CONTROL_CHIP}} \\
  -format BED \\
  -gsize ${{GENOMESIZE}} \\
  -pval {pval} \\
  -shift ${{SHIFT}} \\
  -extsize ${{EXTSIZE}} \\
  -name ${{NAME}} \\
  -output ${{OUTPUT}} \\
  -other "${{OTHER}}" \\
  -dtype "narrowPeak" \\
  -rank ${{RANK}}
echo "Peak-calling done..."
# =========call peak with macs2 (broad)===========
# Parallel running
# =======================
echo "Try to run peak-calling with broad peak mode & sorting peak with {rank}..."
OTHER_BROAD="{otherBroadStr}"
OUTPUT_BROAD="{outputBroadStr}"
callMacs2Parallel.py -cpu ${{THREADS}} \\
  -ip ${{IP_CHIP}} \\
  -input ${{CONTROL_CHIP}} \\
  -format BED \\
  -gsize ${{GENOMESIZE}} \\
  -pval {pval} \\
  -shift ${{SHIFT}} \\
  -extsize ${{EXTSIZE}} \\
  -name ${{NAME}} \\
  -output ${{OUTPUT_BROAD}} \\
  -other "${{OTHER_BROAD}}" \\
  -dtype "broadPeak" \\
  -rank ${{RANK}} \\
echo "Broad-peak-calling done..."
{idrCommand}
'''

idrTemplate = '''
# =========running IDR pipeline ===========
# running on macs2 narrow peaks
# =======================
echo "Running IDR on narrow peaks..."
BLACKLIST="{blacklist}"
IDR_SAMPLE="{idrSampleStr}"
PEAKLIST="{idrPeaklistStr}"
IDR_PREFIX="{idrPrefixStr}"
IDR_OUTPUT="{idrOutputStr}"
IDR_THRESH="{idrThresh}"
RANK="{rank}"
callEncodeIdrParallel.py -cpu ${{THREADS}} \\
  -sample "${{IDR_SAMPLE}}" \\
  -peaklist ${{PEAKLIST}} \\
  -blacklist ${{BLACKLIST}} \\
  -dtype "narrowPeak" \\
  -output ${{IDR_OUTPUT}} \\
  -prefix ${{IDR_PREFIX}} \\
  -thresh ${{IDR_THRESH}} \\
  -rank ${{RANK}}
echo "IDR done (narrow peak)..."
# =========running IDR pipeline ===========
# running on macs2 broad peaks
# =======================
echo "Running IDR on broad peaks..."
IDR_BROAD_SAMPLE="{idrBroadSampleStr}"
BROAD_PEAKLIST="{idrBroadPeaklistStr}"
IDR_BROAD_OUTPUT="{idrBroadOutputStr}"
callEncodeIdrParallel.py -cpu ${{THREADS}} \\
  -sample "${{IDR_BROAD_SAMPLE}}" \\
  -peaklist ${{BROAD_PEAKLIST}} \\
  -blacklist ${{BLACKLIST}} \\
  -dtype "broadPeak" \\
  -output ${{IDR_BROAD_OUTPUT}} \\
  -prefix ${{IDR_PREFIX}} \\
  -thresh ${{IDR_THRESH}} \\
  -rank ${{RANK}}
echo "IDR done (broad peak)..."
'''
# generate sbatch scripts and runBash script
runSbatchScript = os.path.join(bashDir, "runEncodeChipPeakCall.sbatch")
with open(runSbatchScript, 'w') as sbatchO:
    sbatchO.write('#!/bin/sh\n')
    sbatchO.write('BASE={0}\n\n'.format(bashDir))
    sbatchO.write('#sbatch for macs2\n')
    for exp in sorted(expLinkDict.keys()):
        # for macs2: IP, input, name, other output
        # calling peak for narrow peak
        baseExpName = exp
        ipChipStr = ' '.join(expLinkDict[exp]['IP'])
        conChipStr = ' '.join(expLinkDict[exp]['input'])
        nameStr = ' '.join(expLinkDict[exp]['name'])
        otherStr = ','.join(expLinkDict[exp]['other'])
        outputStr = ' '.join(expLinkDict[exp]['output'])
        if bool(args.extsize):
            extsize = args.extsize
        else:
            extsize = int(expInfoDict[exp]['extsize'] / expInfoDict[exp]['repNum'])
        if extsize == 0:
            extsize = 147
        # calling peak for broad peak
        otherBroadStr = ','.join(['--broad --broad-cutoff {0}'.format(pval) 
            for i in range(len(expLinkDict[exp]['name']))])
        outputBroadStr = ' '.join(list(map(lambda x: x.replace('/narrow/', '/broad/'), expLinkDict[exp]['output'])))

        # for idr: input from macs2 narrow peak
        pooledRepName = expLinkDict[exp]['name'][-3]
        pooledRepDir = expLinkDict[exp]['output'][-3]
        repFinalNameList = list(filter(lambda x:re.search(r'final$', x), 
            expLinkDict[exp]['name']))
        repFinalDirList = list(filter(lambda x:re.search(r'final$', x), 
            expLinkDict[exp]['output']))
        repPrNameList = sorted(list(filter(lambda x:re.search(r'rep\d+_pr\d+$', x), 
            expLinkDict[exp]['name'])))
        repPrDirList = sorted(list(filter(lambda x:re.search(r'rep\d+_pr\d+$', x), 
            expLinkDict[exp]['output'])))
        pooledPrNameList = expLinkDict[exp]['name'][-2:]
        pooledPrDirList = expLinkDict[exp]['output'][-2:]
        ## building pairs
        ### true replicates
        idrSampleNameList = [repFinalNameList]
        idrSampleDirList = [repFinalDirList]
        idrPeaklistNameList = [pooledRepName]
        idrPeaklistDirList = [pooledRepDir]
        idrPrefixList = [pooledRepName + '_withRep']
        idrOuputDirList = [pooledRepDir + '_withRep']
        ### self-pseudoreplicates
        repCount = 0
        for i in range(0, len(repPrNameList), 2):
            idrSampleRepPrNameList = [repPrNameList[i], repPrNameList[i+1]]
            idrSampleRepPrDirList = [repPrDirList[i], repPrDirList[i+1]]
            idrSampleNameList.append(idrSampleRepPrNameList)
            idrSampleDirList.append(idrSampleRepPrDirList)
            idrPeaklistNameList.append(repFinalNameList[repCount])
            idrPeaklistDirList.append(repFinalDirList[repCount])
            idrPrefixList.append(repFinalNameList[repCount] + '_withPr')
            idrOuputDirList.append(repFinalDirList[repCount]+ '_withPr')
            repCount += 1
        ### pooled-pseudoreplicates
        idrSampleNameList.append(pooledPrNameList)
        idrSampleDirList.append(pooledPrDirList)
        idrPeaklistNameList.append(pooledRepName)
        idrPeaklistDirList.append(pooledRepDir)
        idrPrefixList.append(pooledRepName + '_withPooledPr')
        idrOuputDirList.append(pooledRepDir + '_withPooledPr')
        ### generate final pairs, narrowPeak
        idrSampleList = list()
        for i in range(len(idrSampleNameList)):
            sample = ' '.join(list(map(lambda x: 
                os.path.join(idrSampleDirList[i][x], idrSampleNameList[i][x] + '_peaks.narrowPeak'), 
                range(len(idrSampleNameList[i])))))
            idrSampleList.append(sample)
        idrPeaklistList = list(map(lambda x: 
            os.path.join(idrPeaklistDirList[x], idrPeaklistNameList[x] + '_peaks.narrowPeak'), 
            range(len(idrPeaklistNameList))))
        idrOuputList = list(map(lambda x:x.replace('/narrow/', '/idr_narrow/'), idrOuputDirList))
        idrSampleStr = ','.join(idrSampleList)
        idrPeaklistStr = ' '.join(idrPeaklistList)
        idrPrefixStr = ' '.join(idrPrefixList)
        idrOutputStr = ' '.join(idrOuputList)
        # for idr: input from macs2 broad peak
        idrBroadSampleList = list(map(lambda x:x.replace('/narrow/', '/broad/'), idrSampleList))
        idrBroadSampleList = list(map(lambda x:x.replace('.narrowPeak', '.broadPeak'), idrBroadSampleList))
        idrBroadPeaklistList = list(map(lambda x:x.replace('/narrow/', '/broad/'), idrPeaklistList))
        idrBroadPeaklistList = list(map(lambda x:x.replace('.narrowPeak', '.broadPeak'), idrBroadPeaklistList))
        idrBroadOuputList = list(map(lambda x:x.replace('/narrow/', '/idr_broad/'), idrOuputDirList))
        idrBroadSampleStr = ','.join(idrBroadSampleList)
        idrBroadPeaklistStr = ' '.join(idrBroadPeaklistList)
        idrBroadOutputStr = ' '.join(idrBroadOuputList)
        # write sbatch file
        logName = "{baseExpName}_PeakCall.macs2_idr.log".format(**vars())
        logPath = os.path.join(logDir, logName)
        sbatchScript = os.path.join(bashDir, baseExpName + '_PeakCall.macs2_idr.sh')
        if args.noIDR:
            idrCommand = ''
        else:
            idrCommand = idrTemplate.format(**vars())
        sbatchCont = sbatchMacs2IdrTemplate.format(**vars())
        with open (sbatchScript, 'w') as sbatchScriptO:
            sbatchScriptO.write(sbatchCont)
        sbatchO.write('sbatch $BASE/{baseExpName}_PeakCall.macs2_idr.sh\n'.format(**vars()))
        os.chmod(sbatchScript, 0o755)
os.chmod(runSbatchScript, 0o755)
