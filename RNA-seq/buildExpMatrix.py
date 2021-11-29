#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input RSEM genes|isoforms counts results directory')
parser.add_argument('-identity', action='store', type=str, choices=['genes', 'isoforms'],
                    default='genes',
                    help='use genes|isoforms to generate exression matrix')
parser.add_argument('-abundance', action='store', type=str, choices=['counts', 'FPKM', 'TPM'],
                    default='counts',
                    help='export types of gene abundance')
parser.add_argument('-grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('-grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output result matrix')
parser.add_argument('-nocqv', action='store_true',
                    default=False,
                    help='output result matrix')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

##public arguments
if bool(args.grepKept):
    kept = True
    regex = re.compile(r'{0}'.format(args.grepKept))
elif bool(args.grepExpel):
    kept = False
    regex = re.compile(r'{0}'.format(args.grepExpel))
else:
    regex = False
countMtxFiles = sorted(glob(os.path.join(args.input, '**', '*.'+args.identity+'.results'), recursive=True))

if args.abundance == 'counts':
    abIndex = 4
elif args.abundance == 'FPKM':
    abIndex = 6
    cqvIndex = 16
else:
    abIndex = 5
    cqvIndex = 13

expDict = defaultdict(list)
cqvDict = defaultdict(list)
sampleNameList = list()
for countMtx in countMtxFiles:
    sampleName = os.path.split(countMtx)[-1].replace('.'+args.identity+'.results', '')
    if bool(regex):
        if kept:
            if bool(regex.search(sampleName)) is False:
                continue
        else:
            if bool(regex.search(sampleName)) is True:
                continue
    sampleNameList.append(sampleName)
    with open(countMtx, 'r') as f:
        __ = f.readline()
        for line in f:
            row = line.strip().split('\t')
            geneId = row[0]
            abundance = row[abIndex]
            expDict[geneId].append(abundance)
            if args.abundance != 'counts':
                cqv = row[cqvIndex]
                cqvDict[geneId].append(cqv)

with open (args.output, 'w') as out:
    if args.identity == 'genes':
        row = ['gene_id']
    else:
        row = ['tx_id']
    row.extend(sampleNameList)
    if args.abundance != 'counts' and args.nocqv is False:
        row.extend(list(map(lambda x: x + '_CQV', sampleNameList)))
    out.write('\t'.join(row) + '\n')
    for geneId in sorted(expDict.keys()):
        row = [geneId]
        row.extend(expDict[geneId])
        if args.abundance != 'counts' and args.nocqv is False:
            row.extend(cqvDict[geneId])
        out.write('\t'.join(row) + '\n')
