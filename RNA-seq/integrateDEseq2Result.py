#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
import math
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input DESeq2 results directory (called by DEseq2Gene.R [*.DESeq2.txt])')
parser.add_argument('-file', nargs='+', type=str,
                    help='input file list instead of -input')
parser.add_argument('-name', nargs='+', type=str,
                    help='designate result names, work with -file')
parser.add_argument('-nindex', nargs='+', type=int,
                    help='extract name as sampleName from filename by index ("_" separated)')
parser.add_argument('-nonLog2', action='store_true',
                    default=False,
                    help='convert log2FC to FC')
parser.add_argument('-grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('-grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output result matrix')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

##public arguments

if bool(args.file) and bool(args.name):
    if len(args.file) != len(args.name):
        parser.print_help()
        parser.exit()

if bool(args.input):
    if bool(args.grepKept):
        kept = True
        regex = re.compile(r'{0}'.format(args.grepKept))
    elif bool(args.grepExpel):
        kept = False
        regex = re.compile(r'{0}'.format(args.grepExpel))
    else:
        regex = False
    deFileList = sorted(glob(os.path.join(args.input, '**', '*.DESeq2.txt'), recursive=True))
elif bool(args.file):
    deFileList = args.file
else:
    parser.print_help()
    parser.exit()

# baseMean
bmDict = defaultdict(str)
# 
deDict = defaultdict(dict)
sampleNameList = list()
for i in range(len(deFileList)):
    deFile = deFileList[i]
    if (bool(args.file) and bool(args.name)) is False:
        sampleName = os.path.split(deFile)[-1]
        if bool(regex):
            if kept:
                if bool(regex.search(sampleName)) is False:
                    continue
            else:
                if bool(regex.search(sampleName)) is True:
                    continue
        sampleName = sampleName.split('.')[0]
        if bool(args.nindex):
            nindexList = args.nindex
            tempRow = sampleName.split('_')
            nameRow = list()
            for nindex in nindexList:
                nameRow.append(tempRow[nindex])
            sampleName = '_'.join(nameRow)
    else:
        sampleName = args.name[i]
    sampleNameList.append(sampleName)
    with open(deFile, 'r') as f:
        __ = f.readline()
        for line in f:
            row = line.strip().split('\t')
            geneId = row[0]
            baseMean = row[1]
            if geneId not in bmDict:
                bmDict[geneId] = baseMean
            fc_val = row[2]
            if (args.nonLog2):
                fc_val = 2 ** float(fc_val)
                fc_val = str(round(12.3456, 4))
            pval = '1' if row[4] == 'NA' else row[4]
            adjp = '1' if row[5] == 'NA' else row[5]
            deDict[geneId][sampleName] = [fc_val, pval, adjp]

with open (args.output, 'w') as out:
    row = ['gene_id', 'baseMean']
    row.extend(['log2FC_' + x for x in sampleNameList])
    row.extend(['pval_' + x for x in sampleNameList])
    row.extend(['adjp_' + x for x in sampleNameList])
    out.write('\t'.join(row) + '\n')
    for geneId in sorted(deDict.keys()):
        row = [geneId, bmDict[geneId]]
        for sampleName in sampleNameList:
            if sampleName not in deDict[geneId]:
                deDict[geneId][sampleName] = ['0', '1', '1']
        row.extend([deDict[geneId][x][0] for x in sampleNameList])
        row.extend([deDict[geneId][x][1] for x in sampleNameList])
        row.extend([deDict[geneId][x][2] for x in sampleNameList])
        out.write('\t'.join(row) + '\n')
