#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input RSEM genes|isoforms expression matrix')
parser.add_argument('-idCol', action='store', type=int,
                    default=0,
                    help='index of gene id column')
parser.add_argument('-nameCol', action='store', type=int,
                    help='index of gene name column')
parser.add_argument('-matrix', action='store', type=str, required=True,
                    help='sample matrix')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output result matrix')
parser.add_argument('-sampleCol', action='store', type=int,
                    default=1,
                    help='index of starting name column of samples')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

colNameList = list()
sampleDict = defaultdict(dict)
with open(args.matrix, 'r') as f:
    colNameList = f.readline().strip().split('\t')[1:]
    for line in f:
        row = line.strip().split('\t')
        sampleName = row[0]
        sampleDict[sampleName] = row[1:]

geneDict = defaultdict(dict)

sampleList = sorted(sampleDict.keys())
expDict = defaultdict(dict)
with open(args.input, 'r') as f:
    allSampleList = f.readline().strip().split('\t')[args.sampleCol:]
    #allSampleList = list(filter(lambda x:(bool(re.search(r'_CQV$', x)) is False), allSampleList))
    expLength = len(allSampleList)
    indexList = list()
    for i in range(expLength):
        if allSampleList[i] in sampleList:
            indexList.append(i)
    for line in f:
        row = line.strip().split('\t')
        geneId = row[args.idCol]
        if bool(args.nameCol):
            geneName = row[args.nameCol]
            geneDict[geneId] = geneName
        expList = row[args.sampleCol:]
        for i in indexList:
            sampleName = allSampleList[i]
            expDict[geneId][sampleName] = expList[i]

with open(args.output, 'w') as out:
    row = ['gene_id', 'sample', 'expression']
    if bool(args.nameCol):
        row = ['gene_id', 'gene_name', 'sample', 'expression']
    row.extend(colNameList)
    out.write('\t'.join(row) + '\n')
    for geneId in sorted(expDict.keys()):
        for sampleName in sorted(expDict[geneId].keys()):
            exp = expDict[geneId][sampleName]
            row = [geneId, sampleName, exp]
            if bool(args.nameCol):
                geneName = geneDict[geneId]
                row = [geneId, geneName, sampleName, exp]
            row.extend(sampleDict[sampleName])
            out.write('\t'.join(row) + '\n')
