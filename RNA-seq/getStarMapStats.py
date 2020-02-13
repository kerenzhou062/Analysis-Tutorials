#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input STAR mapping log file (*.Log.final.out)')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output stats matrix')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

##public arguments
statsDict = defaultdict(int)
with open(args.input, 'r') as f:
    for line in f:
        row = list(map(lambda x:x.strip(), line.strip().split('|')))
        if len(row) < 2:
            continue
        item = row[0]
        value = row[1]
        value = value.replace('%', '')
        statsDict[item] = value

fullCountList = ['Uniquely mapped reads number',
    'Number of reads mapped to multiple loci',
    'Number of reads mapped to too many loci',
    'Number of reads unmapped: too many mismatches',
    'Number of reads unmapped: too short',
    'Number of reads unmapped: other',
    'Number of chimeric reads']

shortList = ['mapped (unique loci)',
    'mapped (multiple loci)',
    'mapped (too many locimapped)',
    'unmapped (too many mismatches)',
    'unmapped (too short)',
    'unmapped (other)',
    'chimeric']

with open(args.output, 'w') as out:
    row = ['type', 'count']
    out.write('\t'.join(row) + '\n')
    for i in range(len(fullCountList)):
        full = fullCountList[i]
        short = shortList[i]
        if full in statsDict:
            value = statsDict[full]
            row = [short, value]
            out.write('\t'.join(row) + '\n')
