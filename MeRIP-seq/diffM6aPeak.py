#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict
import tempfile
import subprocess
import math
from multiprocessing import Pool, Manager

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-c', '--control', action='store', type=str,
                    required=True,
                    help='peak.xls of control group from exomePeak')
parser.add_argument('-t', '--treat', action='store', type=str,
                    required=True,
                    help='peak.xls of treatment group from exomePeak')
parser.add_argument('-d', '--diff', action='store', type=str,
                    help='peak.diff.xls from exomePeak or QNB')
parser.add_argument('--package', action='store', type=str,
                    choices=['exomePeak', 'QNB'],
                    default='QNB',
                    help='package for detection of differential methylated peaks')
parser.add_argument('-p', '--pval', action='store', type=float,
                    default=1,
                    help='cutoff of p-value')
parser.add_argument('-f', '--fold', action='store', type=float,
                    default=1,
                    help='cutoff of fold_enrchment')
parser.add_argument('-q', '--fdr', action='store', type=float,
                    default=1,
                    help='cutoff of fdr')
parser.add_argument('-o', '--output', action='store', type=str,
                    required=True,
                    help='the output peak.custom.diff.xls')
parser.add_argument('-s', '--strand', action='store', type=str,
                    choices=['unstrand', 'reverse', 'forward'],
                    default='reverse',
                    help='library protocols(reverse:fr-firstrand, forward:fr-secondstrand)')
parser.add_argument('--thread', action='store', type=int,
                    help='The number of threads used with --bamdir, run in parallel')
parser.add_argument('--expMtx', action='store', type=str,
                    help='The expression matrix from buildExpMatrix.py')
parser.add_argument('--degMtx', action='store', type=str,
                    help='The differentially expressed gene matrix from DESeq2Gene.R')
parser.add_argument('--grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('--grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('--bamdir', action='store', type=str,
                    help='input sorted bam directory')
parser.add_argument('--cntKey', action='store', type=str,
                    help='keyword for names of control samples')
parser.add_argument('--trtKey', action='store', type=str,
                    help='keyword for names of treatment samples')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def buildBamDict(filter, bamFileList):
    bamDict = defaultdict(list)
    bamRegex = re.compile(r'{0}'.format(filter))
    for bamFile in bamFileList:
        fileName = os.path.split(bamFile)[-1]
        if bool(regex):
            if kept:
                if bool(regex.search(fileName)) is False:
                    continue
            else:
                if bool(regex.search(fileName)) is True:
                    continue
        if bool(bamRegex.search(fileName)) is True:
            if re.search(r'IP', fileName, flags=re.IGNORECASE):
                bamDict['IP'].append(bamFile)
            else:
                bamDict['input'].append(bamFile)
    return bamDict

def runPeakBamRead(peakReadDict, peakFile, bamFile, strand):
    ## get chromosome names from bam.bai
    gsizeTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
    command = 'samtools idxstats {} | cut -f 1 '.format(bamFile)
    chromLineList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
    chromLineList = list(filter(lambda x: (bool(x) and x != '*'), chromLineList))
    with open(gsizeTmp.name, 'w') as temp:
        for line in chromLineList:
            row = line.strip().split('\t')
            temp.write('\t'.join([row[0], '1']) + '\n')
    ## sort peaks by chromosome
    sortPeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
    command = 'bedtools sort -i {} -g {} > {}'.format(peakFile, gsizeTmp.name, sortPeakTmp.name)
    __ = subprocess.check_output(command, shell=True)
    ## intersect peaks with sorted bam
    command = 'bedtools intersect -a {} -b {} -g {} -split -sorted {} -wa -c'.format(sortPeakTmp.name, bamFile, gsizeTmp.name, strand)
    peakReadsList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
    gsizeTmp.close()
    sortPeakTmp.close()
    ## construct peakId-bam-reads
    for line in peakReadsList:
        if bool(line) is False:
            continue
        row = line.strip().split('\t')
        peakId = row[3]
        readsNum = row[-1]
        # the effective way to update nested dict
        add = {bamFile: int(readsNum)}
        peakReadDict[peakId] = dict(peakReadDict[peakId], **add)

##public arguments
if bool(args.grepKept):
    kept = True
    regex = re.compile(r'{0}'.format(args.grepKept))
elif bool(args.grepExpel):
    kept = False
    regex = re.compile(r'{0}'.format(args.grepExpel))
else:
    regex = False

if bool(args.diff):
    if args.package == 'QNB':
        pvalDiff = args.pval
    else:
        pvalDiff = math.log(args.pval, 10)
        fdrDiff = math.log(args.fdr, 10)

pval = math.log(args.pval, 10)
fdr = math.log(args.fdr, 10)
fold = args.fold

controlPeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
treatPeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)

nameRow = list()
with open(args.control, 'r') as f, open(controlPeakTmp.name, 'w') as temp:
    nameRow = f.readline().strip().split('\t') + ["diff.lg.fdr", "diff.lg.p", "diff.log2.fc"]
    for line in f:
        row = line.strip().split('\t')
        if float(row[12]) > pval:
            continue
        if float(row[13]) > fdr:
            continue
        if float(row[14]) < fold:
            continue
        temp.write('\t'.join(row) + '\n')

with open(args.treat, 'r') as f, open(treatPeakTmp.name, 'w') as temp:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        if float(row[12]) > pval:
            continue
        if float(row[13]) > fdr:
            continue
        if float(row[14]) < fold:
            continue
        temp.write('\t'.join(row) + '\n')

command = 'bedtools intersect -a {0} -b {1} -v -s'.format(controlPeakTmp.name, treatPeakTmp.name)
controlUniqPeakList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')

command = 'bedtools intersect -a {0} -b {1} -v -s'.format(treatPeakTmp.name, controlPeakTmp.name)
treatUniqPeakList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
controlPeakTmp.close()
treatPeakTmp.close()

combineRow = list()
for line in controlUniqPeakList:
    if bool(line) is False:
        continue
    row = line.strip().split('\t')
    pvalLog = row[12]
    fdrLog = row[13]
    log2fc = "-99"
    row = row + [fdrLog, pvalLog, log2fc]
    combineRow.append(row)

for line in treatUniqPeakList:
    if bool(line) is False:
        continue
    row = line.strip().split('\t')
    pvalLog = row[12]
    fdrLog = row[13]
    log2fc = "99"
    row = row + [fdrLog, pvalLog, log2fc]
    combineRow.append(row)

if bool(args.diff):
    combinePeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
    diffPeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
    with open(combinePeakTmp.name, 'w') as temp:
        for row in combineRow:
            temp.write('\t'.join(row) + '\n')
    with open(args.diff, 'r') as f, open(diffPeakTmp.name, 'w') as temp:
        __ = f.readline()
        for line in f:
            temp.write(line)
    
    command = 'bedtools intersect -a {0} -b {1} -v -s'.format(diffPeakTmp.name, combinePeakTmp.name)
    diffUniqPeakList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
    combinePeakTmp.close()
    diffPeakTmp.close()
    
    for line in diffUniqPeakList:
        if bool(line) is False:
            continue
        row = line.strip().split('\t')
        if args.package == 'QNB':
            if row[15] == 'NA':
                continue
            elif float(row[15]) <= pvalDiff:
                row[15] = str(math.log(float(row[15]), 10))
                row[17] = row[16]
                row[16] = row[15]
                combineRow.append(row)
        else:
            if float(row[16]) <= pvalDiff and float(row[15]) <= fdrDiff:
                combineRow.append(row)

# build expression dict
expSampleDict = defaultdict(list)
expDict = defaultdict(dict)
if bool(args.expMtx):
    with open(args.expMtx, 'r') as f:
        expNameRow = f.readline().split('\t')
        rowLength = len(expNameRow)
        for i in range(1,rowLength):
            if bool(args.cntKey) and bool(args.trtKey):
                if bool(re.search(r'{}'.format(args.cntKey), expNameRow[i])):
                    expSampleDict['cnt'].append(i)
                if bool(re.search(r'{}'.format(args.trtKey), expNameRow[i])):
                    expSampleDict['trt'].append(i)
            else:
                sys.stderr.write('--cntKey and --trtKey required by --expMtx!')
                sys.exit()
        for line in f:
            row = line.strip().split('\t')
            geneId = row[0].split('.')[0]
            expDict[geneId]['cnt'] = sum(map(lambda x: float(row[x]), expSampleDict['cnt'])) / len(expSampleDict['cnt'])
            expDict[geneId]['trt'] = sum(map(lambda x: float(row[x]), expSampleDict['trt'])) / len(expSampleDict['trt'])

# build differentlly expressed gene dict
deGeneDict = defaultdict(dict)
if bool(args.degMtx):
    with open(args.degMtx, 'r') as f:
        __ = f.readline()
        for line in f:
            row = line.strip().split('\t')
            geneId = row[0].split('.')[0]
            log2fc = '{0:.2f}'.format(float(row[2]))
            padj = row[-1]
            if padj == 'NA':
                padj = '1'
            else:
                padj = '{0:.2e}'.format(float(row[-1]))
            deGeneDict[geneId] = [log2fc, padj]

peakIdList = list()
peakDict = defaultdict(dict)
for i in range(len(combineRow)):
    geneId = combineRow[i][3].split('.')[0]
    combineRow[i][3] = '|'.join([combineRow[i][3], str(i)])
    peakId = combineRow[i][3]
    peakDict[peakId]['row'] = combineRow[i]
    peakIdList.append(peakId)
    ## initiate bam reads
    peakDict[peakId]['reads'] = list()
    ## add expression
    if bool(expDict):
        expList = ['NA', 'NA']
        if 'cnt' in expDict[geneId]:
            expList[0] = '{0:.2f}'.format(expDict[geneId]['cnt'])
        if 'trt' in expDict[geneId]:
            expList[1] = '{0:.2f}'.format(expDict[geneId]['trt'])
        peakDict[peakId]['exp'] = expList
    else:
        peakDict[peakId]['degene'] = list()
    ## add differential expression
    if bool(deGeneDict):
        if geneId not in deGeneDict:
            peakDict[peakId]['degene'] = ['0', '1']
        else:
            peakDict[peakId]['degene'] = deGeneDict[geneId]
    else:
        peakDict[peakId]['degene'] = list()

if bool(args.bamdir):
    ## calculate reads from samples that cover peaks
    ## innitiate variables
    if bool(args.thread) is False:
        try:
            thread = int(os.environ['SLURM_JOB_CPUS_PER_NODE'].split('(')[0]) * int(os.environ['SLURM_JOB_NUM_NODES'])
        except KeyError as e:
            thread = os.cpu_count()
    else:
        thread = args.thread
    ## run in parallel
    pool = Pool(processes=thread)
    peakReadDict = Manager().dict()
    cntBamDict = defaultdict(list)
    trtBamDict = defaultdict(list)
    if args.strand == 'unstrand':
        strand = ''
    elif args.strand == 'reverse':
        strand = '-S'
    else:
        strand = '-s'
    bamFileList = sorted(glob(os.path.join(args.bamdir, '**', '*.bam'), recursive=True))
    if bool(args.cntKey) and bool(args.trtKey):
        cntBamDict = buildBamDict(args.cntKey, bamFileList)
        trtBamDict = buildBamDict(args.trtKey, bamFileList)
    else:
        sys.stderr.write('--cntKey and --trtKey required by --bamdir!')
        sys.exit()

    peakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
    with open(peakTmp.name, 'w') as temp:
        for row in combineRow:
            peakId = row[3]
            peakReadDict[peakId] = dict()
            temp.write('\t'.join(row) + '\n')

    ## run peaks-bam-reads in parallel
    bamFileList = cntBamDict['IP'] + cntBamDict['input']  + trtBamDict['IP'] + trtBamDict['input']
    for bamFile in bamFileList:
        #runPeakBamRead(peakReadDict, peakTmp.name, bamFile, strand)
        pool.apply_async(runPeakBamRead, args=(peakReadDict, peakTmp.name, bamFile, strand))
    pool.close()
    pool.join()
    peakTmp.close()
    ## calculate real.diff.log2.fc and average reads of samples
    bamList = [cntBamDict['IP'], cntBamDict['input'], trtBamDict['IP'], trtBamDict['input']]
    for peakId in peakIdList:
        valueList = ['NA' for i in range(4)]
        for i in range(4):
            if bool(bamList[i]):
                valueList[i] = sum(map(lambda x: peakReadDict[peakId][x], bamList[i])) / len(bamList[i])
        if 'NA' in valueList:
            realLog2FC = 'NA'
        else:
            realFC = (valueList[2] + 1) /(valueList[3] + 1) * (valueList[1] + 1) / (valueList[0] + 1)
            realLog2FC = math.log(realFC, 2)
        valueRow = [realLog2FC] + valueList
        valueRow = list(map(lambda x: str('{0:.2f}'.format(x)), valueRow))
        peakDict[peakId]['reads'] = valueRow

with open(args.output, 'w') as out:
    addNameRow = list()
    if bool(args.bamdir):
        nameRow.extend(['real.diff.log2.fc', 'cntIpAve', 'cntInputAve', 'trtIpAve', 'trtInputAve'])
    if bool(args.expMtx):
        nameRow.extend(['cntExpAve', 'trtExpAve'])
    if bool(args.degMtx):
        nameRow.extend(['DE.log2fc', 'DE.fdr'])
    out.write('\t'.join(nameRow) + '\n')
    for peakId in peakIdList:
        row = peakDict[peakId]['row']
        row[3] = row[3].split('|')[0]
        row += peakDict[peakId]['reads'] + peakDict[peakId]['exp'] + peakDict[peakId]['degene']
        out.write('\t'.join(row) + '\n')
