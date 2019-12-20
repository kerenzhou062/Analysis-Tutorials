#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict

#usage: runExomePeakBash.py or runExomePeakBash.py <bam dir>

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input RSEM genes|isoforms counts results directory')
parser.add_argument('-identity', action='store', type=str, choices=['genes', 'isoforms'],
                    default='genes',
                    help='use genes|isoforms to generate exression matrix')
parser.add_argument('-grep', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output result matrix')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

##public arguments
if bool(args.grep):
    regex = re.compile(r'{0}'.format(args.grep))
countMtxFiles = sorted(glob(os.path.join(args.input, '**', '*.'+args.identity+'.results'), recursive=True))

expDict = defaultdict(list)
sampleNameList = list()
for countMtx in countMtxFiles:
    sampleName = os.path.split(countMtx)[-1].replace('.'+args.identity+'.results', '')
    if bool(args.grep):
        if bool(regex.search(sampleName)) is False:
            continue
    sampleNameList.append(sampleName)
    with open(countMtx, 'r') as f:
        __ = f.readline()
        for line in f:
            row = line.strip().split('\t')
            geneId = row[0]
            exp_counts = row[4]
            expDict[geneId].append(exp_counts)

with open (args.output, 'w') as out:
    row = ['gene_id']
    row.extend(sampleNameList)
    out.write('\t'.join(row) + '\n')
    for geneId in sorted(expDict.keys()):
        row = [geneId]
        row.extend(expDict[geneId])
        out.write('\t'.join(row) + '\n')
