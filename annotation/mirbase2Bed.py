#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict

#usage: mirBase2Bed.py -inpupt miRNA.v22.gff3 --tx -output miRNA.v22.bed12

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str,
                     help='input GFF3 file')
parser.add_argument('--tx', action='store_true',
                    default=False, help='turn on transcript isoform mode')
parser.add_argument('-output', action='store', type=str,
                    help='output bed')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

with open(args.input, 'r') as f, open(args.output, 'w') as out:
    for line in f:
        if re.match(r'^#', line):
            continue
        row = line.strip().split('\t')
        feature = row[2]
        chrom = row[0]
        start = str(int(row[3]) - 1)
        end = row[4]
        strand = row[6]
        infoDict = defaultdict(dict)
        for item in row[8].split(';'):
            key = item.split('=')[0]
            value = item.split('=')[1]
            infoDict[key] = value
        # output
        col4th = ''
        if args.tx:
            if feature == 'miRNA':
                geneId = infoDict['Derives_from']
                txID = infoDict['ID']
                txName = infoDict['Name']
                col4th = ':'.join([txID, txName, geneId, feature])
            else:
                priTxID = infoDict['ID']
                priTxName = infoDict['Name']
                col4th = ':'.join([priTxID, priTxName, 'NA', feature])
        else:
            if feature == 'miRNA':
                continue
            priTxID = infoDict['ID']
            priTxName = infoDict['Name']
            col4th = ':'.join([priTxID, priTxName, 'NA', feature])
        row = [chrom, start, end, col4th, '0', strand]
        out.write('\t'.join(row) + '\n')
