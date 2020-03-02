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
parser.add_argument('--grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('--grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('--bamdir', action='store', type=str,
                    help='input sorted bam directory')
parser.add_argument('--cntbam', action='store', type=str,
                    help='keyword for bams of control samples')
parser.add_argument('--trtbam', action='store', type=str,
                    help='keyword for bams of treatment samples')

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

if bool(args.bamdir):
    cntBamDict = defaultdict(list)
    trtBamDict = defaultdict(list)
    bamFileList = sorted(glob(os.path.join(args.bamdir, '**', '*.bam'), recursive=True))
    if bool(args.cntbam):
        cntBamDict = buildBamDict(args.cntbam, bamFileList)
    if bool(args.trtbam):
        trtBamDict = buildBamDict(args.trtbam, bamFileList)

    peakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
    peakIdList = list()
    peakReadsDict = defaultdict(dict)
    with open(peakTmp.name, 'w') as temp:
        for i in range(len(combineRow)):
            row = combineRow[i]
            row[3] = '|'.join([row[3], str(i)])
            peakReadsDict[row[3]]['row'] = row
            peakReadsDict[row[3]]['reads'] = defaultdict(int)
            peakIdList.append(row[3])
            temp.write('\t'.join(row) + '\n')

    bamFileList = cntBamDict['IP'] + cntBamDict['input']  + trtBamDict['IP'] + trtBamDict['input']
    for i in range(len(bamFileList)):
        bamFile = bamFileList[i]
        gsizeTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
        command = 'samtools idxstats {} | cut -f 1 '.format(bamFile)
        chromLineList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
        chromLineList = list(filter(lambda x: (bool(x) and x != '*'), chromLineList))
        with open(gsizeTmp.name, 'w') as temp:
            for line in chromLineList:
                row = line.strip().split('\t')
                temp.write('\t'.join([row[0], '1']) + '\n')
        sortPeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
        command = 'bedtools sort -i {} -g {} > {}'.format(peakTmp.name, gsizeTmp.name, sortPeakTmp.name)
        __ = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
        command = 'bedtools intersect -a {} -b {} -g {} -split -sorted -s -wa -c'.format(sortPeakTmp.name, bamFile, gsizeTmp.name)
        peakReadsList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
        peakReadsList = list(filter(lambda x:bool(x), peakReadsList))
        for line in peakReadsList:
            row = line.strip().split('\t')
            peakId = row[3]
            readsNum = row[-1]
            peakReadsDict[peakId]['reads'][bamFile] = int(readsNum)
        gsizeTmp.close()
        sortPeakTmp.close()
    peakTmp.close()
    bamList = [cntBamDict['IP'], cntBamDict['input'], trtBamDict['IP'], trtBamDict['input']]
    nameList = ['real.diff.log2.fc', 'cntIpAve', 'cntInputAve', 'trtIpAve', 'trtInputAve']
    with open(args.output, 'w') as out:
        nameRow = nameRow + nameList
        out.write('\t'.join(nameRow) + '\n')
        for peakId in peakIdList:
            peakRow = peakReadsDict[peakId]['row']
            peakRow[3] = peakRow[3].split('|')[0]
            valueList = ['NA' for i in range(4)]
            tempDict = peakReadsDict[peakId]['reads']
            for i in range(len(bamList)):
                if bool(bamList[i]):
                    valueList[i] = sum(map(lambda x: tempDict[x], bamList[i])) / len(bamList[i])
            if 'NA' in valueList:
                realFC = 'NA'
            else:
                realFC = (valueList[2] + 1) /(valueList[3] + 1) * (valueList[1] + 1) / (valueList[0] + 1)
                realLog2FC = math.log(realFC, 2)
            valueRow = [realLog2FC] + valueList
            valueRow = list(map(lambda x: str('{0:.2f}'.format(x)), valueRow))
            row = peakRow + valueRow
            out.write('\t'.join(row) + '\n')
else:
    with open(args.output, 'w') as out:
        out.write('\t'.join(nameRow) + '\n')
        for row in combineRow:
            out.write('\t'.join(row) + '\n')
