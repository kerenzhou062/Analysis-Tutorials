#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str, required=True,
                    help='results of intersection between miCLIP bed and m6A bedgraph')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output result name')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

nullId = list()
sampleDict = defaultdict(int)
siteDict = defaultdict(dict)
with open(args.input, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        siteId = row[3]
        bgFileName = row[12]
        if bgFileName == '.':
            nullId.append(siteId)
        else:
            sample = row[12].split('.')[0]
            expReadCount = row[17]
            siteDict[siteId][sample] = str(round(float(expReadCount), 0))
            sampleDict[sample] += 1

sampleList = sorted(sampleDict.keys())
with open(args.output, 'w') as out:
    row = ['gene_id']
    row.extend(sampleList)
    out.write('\t'.join(row) + '\n')
    for siteId in nullId:
        row = [siteId]
        row.extend(['0' for x in sampleList])
        out.write('\t'.join(row) + '\n')
    for siteId in sorted(siteDict.keys()):
        row = [siteId]
        for sample in sampleList:
            if sample in siteDict[siteId]:
                row.append(siteDict[siteId][sample])
            else:
                row.append('0')
        out.write('\t'.join(row) + '\n')
