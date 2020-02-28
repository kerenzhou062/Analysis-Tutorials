#!/usr/bin/env python3
import os
import sys
import argparse
import re
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

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

##public arguments
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

combineRow = [nameRow]
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

with open(args.output, 'w') as out:
    for row in combineRow:
        out.write('\t'.join(row) + '\n')
