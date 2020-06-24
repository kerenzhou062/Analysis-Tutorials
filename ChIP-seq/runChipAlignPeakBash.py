#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict

# usage: runChipSeBowtie2AlignBash.py
# or runChipSeBowtie2AlignBash.py <fastq dir>

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-alignment', action='store', type=str,
                    default='alignment',
                    help='alignment result directory')
parser.add_argument('-basepath', action='store', type=str,
                    default='./',
                    help='base directory for -log, -bash, -alignment, -peak')
parser.add_argument('-bash', action='store', type=str,
                    default='runBash',
                    help='runBash directory')
parser.add_argument('-blacklist', action='store', type=str,
                    required=True,
                    help='blacklist for genome: (prefix or file)')
parser.add_argument('-cpu', action='store', type=int,
                    default=10,
                    help='threads used for bowtie2')
parser.add_argument('-chrsize', action='store', type=str,
                    default=10,
                    help='chromosome size file (eg. hg38.chrom.sizes)')
parser.add_argument('-extsize', action='store', type=int,
                    help='--extsize parameter (macs2, histone:147)')
parser.add_argument('-fasta', action='store', type=str,
                    required=True, 
                    help='fasta for genome')
parser.add_argument('-gsize', action='store', type=str,
                    required=True,
                    help='-g parameter (macs2, mm|hs|dm|ce, or exact number)')
parser.add_argument('-idrcutoff', action='store', type=float,
                    default=0.05, help='--soft-idr-threshold parameter (idr)')
parser.add_argument('-index', action='store', type=str,
                    required=True,
                    help='bowtie2 genome index (prefix or file)')
parser.add_argument('-input', action='store', type=str,
                    required=True,
                    help='input fastq directory')
parser.add_argument('-log', action='store', type=str,
                    default='log',
                    help='log directory of runBash')
parser.add_argument('-matrix', action='store', type=str,
                    required=True,
                    help='input matrix file (./workflow/sample.matrix)')
parser.add_argument('-memory', action='store', type=str,
                    default='100G',
                    help='memory used for sbatch')
parser.add_argument('-noIDR', action='store_true',
                    default=False,
                    help='disable idr calling')
parser.add_argument('-part', action='store', type=str,
                    default='all',
                    help='partition of slurm server')
parser.add_argument('-peak', action='store', type=str,
                    default='peak',
                    help='peak results directory')
parser.add_argument('-pval', action='store', type=str,
                    default='1e-2', help='-p parameter (macs2)')
parser.add_argument('-quality', action='store', type=int,
                    default=30,
                    help='the mapping quality for "plotFingerprint"')
parser.add_argument('-rank', action='store', type=str,
                    choices=['p.value','q.value','signal.value'],
                    default='p.value', help='--rank parameter (idr)')
parser.add_argument('-shift', action='store', type=str,
                    default='0', help='--shift parameter (macs2)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def CreateDir(root, directory):
    dirPath = os.path.join(root, directory)
    dirPath = os.path.realpath(dirPath)
    os.makedirs(dirPath, exist_ok=True)
    return dirPath

# public arguments
thread = args.cpu
memory = args.memory
partition = args.part
quality = args.quality
gsize = args.gsize
idrcutoff = args.idrcutoff
memory = args.memory
partition = args.part
pval = args.pval
rank = args.rank
shift = args.shift
blacklist = args.blacklist
chrsize = args.chrsize
fastqDir = args.input

## required directory and file
basepath = os.path.realpath(args.basepath)
genomeIndex = args.index
blacklist = args.blacklist
fasta = args.fasta
index = args.index

## preset arguments
genome = os.path.basename(args.index)

bashDir = CreateDir(basepath, args.bash)
mainAlignDir = CreateDir(basepath, args.alignment)
mainPeakDir = CreateDir(basepath, args.peak)
logDir = CreateDir(basepath, args.log)

# build exp->file relationship dictionary
expFqDict = defaultdict(dict)
with open(args.matrix, 'r') as f:
    ## prefix cell expAcc target runType rep1,rep2 IP_rep1_run1:IP_rep1_run2,IP_rep2 input_rep1_run1:input_rep1_run2,input_rep2
    ## IP_rep1_run1 = IP_rep1_run1_R1;IP_rep1_run1_R2
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        exp = row[0]
        runType = row[4]
        repList = row[5].split(',')
        ipFileList = row[6].split(',')
        inputFileList = row[7].split(',')
        ipNameList = ['IP', 'input']
        fileList = [ipFileList, inputFileList]
        for i in range(len(ipNameList)):
            ip = ipNameList[i]
            for j in range(len(repList)):
                rep = repList[j]
                runFileList = fileList[i][j].split(':')
                for runFile in runFileList:
                    pairedFileList = runFile.split(';')
                    for k in range(len(pairedFileList)):
                        pairedNum = k + 1
                        basename = pairedFileList[k]
                        if ip in expFqDict[exp]:
                            if pairedNum in expFqDict[exp][ip][rep]:
                                expFqDict[exp][ip][rep][pairedNum].append(basename)
                            else:
                                expFqDict[exp][ip][rep][pairedNum] = [basename]
                        else:
                            expFqDict[exp][ip] = defaultdict(dict)
                            expFqDict[exp][ip][rep][pairedNum] = [basename]
# ================================
# sbatch main template:
# step 1 to 2c, 2g
# ================================
chipAlignPeakCallTemplate = '''#!/bin/bash
#SBATCH --job-name={exp}_ChipAlign_    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Mail user
#SBATCH -n {thread}                # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p {partition}                # default queue is all if you don't specify
#SBATCH --mem={memory}                # Amount of memory in GB
#SBATCH --time=120:00:00               # Time limit hrs:min:sec
#SBATCH --output={logFile}   # Standard output and error log

### this pipeline mainly based on ENCODE3 pipeline:
###https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c

# ================================
# load required module
# ================================
module load picard/2.21.1
module load bedtools/2.29.0
module load Python/3.6.5
module load Anaconda2/4.2.0

# ===========step 1==============
# Align fastqs to genome
# ================================
echo "align fastqs to genome..."
FASTQ_DIR={fastqDir}

{alignCommand}

echo "Reads alignment done!"

# get or set extsize
echo "get or set extsize"
{getExtsizeCommand}

echo "Start to call peaks & run IDR"
{peakCallCommand}

echo "Start to overlap peaks"
{peakOverlapCommand}

echo "Pipeline done!"
'''

runAlignSingleTemplate='''
echo "Algin fastq of {alignPrefix} to genome:{genome}"

if [ ! -d "{algnRoot}" ]; then
  mkdir -p "{algnRoot}"
fi

alignLog="{alignLog}"
runChipAlignGenome.sh --prefix {alignPrefix} --quality {quality} \\
  --thread {thread} --trim 50 \\
  --fasta {fasta} \\
  --index {index} \\
  --root {algnRoot} \\
  --output {alignOutput} \\
  --fastq1 {fastq1} \\
  > $alignLog 2>&1
'''

runAlignPairedTemplate='''
echo "Algin fastq of {alignPrefix} to genome:{genome}"

if [ ! -d "{algnRoot}" ]; then
  mkdir -p "{algnRoot}"
fi

alignLog="{alignLog}"
runChipAlignGenome.sh --prefix {alignPrefix} --quality {quality} \\
  --thread {thread} --trim 50 \\
  --fasta {fasta} \\
  --index {index} \\
  --root {algnRoot} \\
  --output {alignOutput} \\
  --fastq1 {fastq1} \\
  --fastq2 {fastq2} \\
  > $alignLog 2>&1
'''

runChipPoolReplicatesTemplate='''
echo "Run runChipPoolReplicates for {exp} ({ip})"
chipPoolRepLog="{chipPoolRepLog}"
runChipPoolReplicates.sh --prefix {ipPrefix} \\
  --blacklist {blacklist} \\
  --input {algnRoot} \\
  > $chipPoolRepLog 2>&1
'''

runChipPlotJsdTemplate='''
echo "Run runChipPlotJsd for {exp}"
chipPlotJdLog="{chipPlotJdLog}"
runChipPlotJsd.sh --thread {thread} --quality {quality} \\
  --prefix {exp} \\
  --input {expDir} \\
  > $chipPlotJdLog 2>&1
'''

macs2PeakTemplate = '''
BASE="{basepath}"
MAIN_ALIGN_DIR="{mainAlignDir}"
MAIN_PEAK_DIR="{mainPeakDir}"
#separated by " "
THREADS={thread}
GENOMESIZE="{gsize}"
IP_CHIP="{ipChipStr}"
CONTROL_CHIP="{conChipStr}"
NAME="{nameStr}"
OTHER="{otherStr}"
OUTPUT="{outputStr}"
SHIFT="{shift}"
CHRSIZE="{chrsize}"
RANK="{rank}"
echo "Change wd directory to ${{MAIN_PEAK_DIR}}"
cd ${{MAIN_PEAK_DIR}}

# =========call peak with macs2 ===========
# Parallel running
# =======================
echo "Peak-calling starting & sorting narrowPeak with {rank}..."

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

echo "narrowPeak-calling done..."

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

echo "broadPeak-calling done..."
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
IDR_THRESH="{idrcutoff}"
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

peakOverlapTemplate = '''
# =========overlap peaks on real pealicates===========
echo "Running overlap peaks on {peakType} peaks..."
runChipOverlapPeaks.sh --type {peakType} --prefix {exp} \\
  --output {overlapDir} \\
  --blacklist {blacklist} \\
  --peak1 {peak1} \\
  --peak2 {peak2}
echo "overlap {peakType} done..."
'''

# generate sbatch scripts and runBash script
runSbatchScript = os.path.join(bashDir, "runChipPipleLine.sbatch")
with open(runSbatchScript, 'w') as sbatch:
    sbatch.write('#!/bin/sh\n')
    sbatch.write('#sbatch for ChIP-seq samples\n')
    sbatch.write('BASE=' + bashDir + '\n\n')
    # expFqDict: exp->IP->rep->[fq1,fq2]
    for exp in sorted(expFqDict.keys()):
        runAlignCommandList = list()
        ipList = sorted(expFqDict[exp].keys())
        if len(ipList) == 1:
            sys.error.write("Lack of IP or input for " + exp + '\n')
            sys.exit()
        for ip in ipList:
            ipDict = expFqDict[exp][ip]
            algnRoot = os.path.join(mainAlignDir, exp, ip)
            ipPrefix = '_'.join([exp, ip])
            for rep in sorted(ipDict.keys()):
                repDict = ipDict[rep]
                pairedNumList = sorted(repDict.keys())
                pairedNum = len(pairedNumList)
                ## bash arguments
                alignPrefix = '_'.join([exp, ip, rep])
                alignLog = os.path.join(algnRoot, alignPrefix + '.log')
                alignOutput = os.path.join(algnRoot, rep)
                fastq1 = ','.join(map(lambda x: os.path.join("$FASTQ_DIR/", x) , repDict[1]))
                ## runChipAlignGenome command
                if pairedNum == 1:
                    runAlignCommandList.append(runAlignSingleTemplate.format(**vars()))
                else:
                    fastq2 = ','.join(map(lambda x: os.path.join("$FASTQ_DIR/", x) , repDict[2]))
                    runAlignCommandList.append(runAlignPairedTemplate.format(**vars()))
            ## runChipPoolReplicates command
            chipPoolRepLog = os.path.join(algnRoot, ipPrefix + '.ChipPoolReplicates.log' )
            runAlignCommandList.append(runChipPoolReplicatesTemplate.format(**vars()))
        ## ## runChipPoolReplicates command
        expDir = os.path.join(mainAlignDir, exp)
        chipPlotJdLog = os.path.join(expDir, exp + '.ChipPlotJd.log' )
        runAlignCommandList.append(runChipPlotJsdTemplate.format(**vars()))
        alignCommand = '\n'.join(runAlignCommandList)
        ## get extsize from *.trim.filt.sample.15.tagAlign.gz.cc.qc of IP samples
        expRoot = os.path.join(mainAlignDir, exp)
        if bool(args.extsize) is True:
            extsize = args.extsize
            getExtsizeCommand = 'EXTSIZE={extsize}'.format(**vars())
        else:
            getExtsizeCommand = 'EXTSIZE=$(getEstFragLen.py -input {expRoot})'.format(**vars())
        ## call peaks command
        tagAlignFinalApp = '.filt.srt.nodup.final.tagAlign.gz'
        tagAlignPr1App = '.filt.nodup.pr1.tagAlign.gz'
        tagAlignPr2App = '.filt.nodup.pr2.tagAlign.gz'
        tagAlignPoolPr1App = '.pooled.pr1.tagAlign.gz'
        tagAlignPoolPr2App = '.pooled.pr2.tagAlign.gz'
        tagAlignPoolApp = '.pooled.tagAlign.gz'
        ## get information from IP
        expLinkDict = defaultdict(list)
        repList = sorted(expFqDict[exp]['IP'])
        inputRepList = sorted(expFqDict[exp]['input'])
        repNum = len(repList)
        if repNum != len(inputRepList):
            sys.error.write("Not equal IP and input for " + exp + '\n')
            continue
        if repNum > 2:
            sys.error.write("This script does not support more than 2 relicates: " + exp + '\n')
            continue
        expLinkDict['name'].extend(list(map(lambda x:'_'.join([exp, x, 'final']), repList)))
        expLinkDict['name'].extend(list(map(lambda x:'_'.join([exp, x, 'pr1']), repList)))
        expLinkDict['name'].extend(list(map(lambda x:'_'.join([exp, x, 'pr2']), repList)))
        expLinkDict['name'].append('_'.join([exp, 'pooled']))
        expLinkDict['name'].append('_'.join([exp, 'pr1', 'pooled']))
        expLinkDict['name'].append('_'.join([exp, 'pr2', 'pooled']))
        expLinkDict['other'] = ["--keep-dup all -B --SPMR" for i in range(len(expLinkDict['name'])) ]
        outputList = list(map(lambda x:os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', \
            '_'.join([exp, x, 'final'])), repList))
        outputList.extend(list(map(lambda x:os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', \
            '_'.join([exp, x, 'pr1'])), repList)))
        outputList.extend(list(map(lambda x:os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', \
            '_'.join([exp, x, 'pr2'])), repList)))
        outputList.append(os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', exp+'_pooled'))
        outputList.append(os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', exp+'_pr1_pooled'))
        outputList.append(os.path.join( "${MAIN_PEAK_DIR}", exp, 'narrow', exp+'_pr2_pooled'))
        expLinkDict['output'] = outputList
        for ip in ipList:
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
            expLinkDict[ip] = ipTagAlignList
        ## construct peak-call arguments
        ipChipStr = ' '.join(expLinkDict['IP'])
        conChipStr = ' '.join(expLinkDict['input'])
        nameStr = ' '.join(expLinkDict['name'])
        otherStr = ','.join(expLinkDict['other'])
        outputStr = ' '.join(expLinkDict['output'])
        ### calling peak for broad peak
        otherBroadStr = ','.join(['--broad --broad-cutoff {0}'.format(pval) 
            for i in range(len(expLinkDict['name']))])
        outputBroadStr = ' '.join(list(map(lambda x: x.replace('/narrow/', '/broad/'), expLinkDict['output'])))
        ### for idr: input from macs2 narrow peak
        pooledRepName = expLinkDict['name'][-3]
        pooledRepDir = expLinkDict['output'][-3]
        repFinalNameList = list(filter(lambda x:re.search(r'final$', x), 
            expLinkDict['name']))
        repFinalDirList = list(filter(lambda x:re.search(r'final$', x), 
            expLinkDict['output']))
        repPrNameList = sorted(list(filter(lambda x:re.search(r'rep\d+_pr\d+$', x), 
            expLinkDict['name'])))
        repPrDirList = sorted(list(filter(lambda x:re.search(r'rep\d+_pr\d+$', x), 
            expLinkDict['output'])))
        pooledPrNameList = expLinkDict['name'][-2:]
        pooledPrDirList = expLinkDict['output'][-2:]
        ### building pairs
        #### true replicates
        idrSampleNameList = [repFinalNameList]
        idrSampleDirList = [repFinalDirList]
        idrPeaklistNameList = [pooledRepName]
        idrPeaklistDirList = [pooledRepDir]
        idrPrefixList = [pooledRepName + '_withRep']
        idrOuputDirList = [pooledRepDir + '_withRep']
        #### self-pseudoreplicates
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
        #### pooled-pseudoreplicates
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
        ### for idr: input from macs2 broad peak
        idrBroadSampleList = list(map(lambda x:x.replace('/narrow/', '/broad/'), idrSampleList))
        idrBroadSampleList = list(map(lambda x:x.replace('.narrowPeak', '.broadPeak'), idrBroadSampleList))
        idrBroadPeaklistList = list(map(lambda x:x.replace('/narrow/', '/broad/'), idrPeaklistList))
        idrBroadPeaklistList = list(map(lambda x:x.replace('.narrowPeak', '.broadPeak'), idrBroadPeaklistList))
        idrBroadOuputList = list(map(lambda x:x.replace('/narrow/', '/idr_broad/'), idrOuputDirList))
        idrBroadSampleStr = ','.join(idrBroadSampleList)
        idrBroadPeaklistStr = ' '.join(idrBroadPeaklistList)
        idrBroadOutputStr = ' '.join(idrBroadOuputList)
        if args.noIDR:
            idrCommand = ''
        else:
            idrCommand = idrTemplate.format(**vars())
        peakCallCommand = macs2PeakTemplate.format(**vars())
        ## overlap peaks on real replicates
        peakOverlapCommand = ''
        if repNum == 2:
            peakOverlapCommandList = list()
            ### eg: Kas1_control_NA_rep1_final/Kas1_control_NA_rep1_final_peaks.narrowPeak
            for peakType in ['narrowPeak', 'broadPeak']:
                overlapDir = os.path.join(mainPeakDir, exp, 'overlap')
                if peakType == 'narrowPeak':
                    peakDir = os.path.join(mainPeakDir, exp, 'narrow')
                else:
                    peakDir = os.path.join(mainPeakDir, exp, 'broad')
                peak1DirName = '_'.join([exp, repList[0], 'final'])
                peak2DirName = '_'.join([exp, repList[1], 'final'])
                peak1Name = '_'.join([exp, repList[0], 'final_peaks.']) + peakType
                peak2Name = '_'.join([exp, repList[1], 'final_peaks.']) + peakType
                peak1 = os.path.join(peakDir, peak1DirName, peak1Name)
                peak2 = os.path.join(peakDir, peak2DirName, peak2Name)
                peakOverlapCommandList.append(peakOverlapTemplate.format(**vars()))
            peakOverlapCommand = '\n'.join(peakOverlapCommandList)
        ## final script
        logFile = os.path.join(logDir, 'ChipAlignPeakCall.' + exp + '.log')
        expBashFileName = 'ChipAlignPeakCall.' + exp + '.sh'
        expBashFile = os.path.join(bashDir, expBashFileName)
        expChipAlignPeakCallMainBash = chipAlignPeakCallTemplate.format(**vars())
        with open(expBashFile, 'w') as expBash:
            expBash.write(expChipAlignPeakCallMainBash)
        sbatch.write('sbatch $BASE/{expBashFileName}\n'.format(**vars()))
        os.chmod(expBashFile, 0o755)
os.chmod(runSbatchScript, 0o755)
