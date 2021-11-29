#!/usr/bin/env python3
import os
import sys
import argparse
import re


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input fasta file')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output result (coordiantes in 0-base)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

chromRegex = re.compile(r'^>')
chromKeepRegex = re.compile(r'^>chr\d+|>chrX|>chrY|>chrM')
chrFaDict = dict()

chrom = ''
with open(args.input, 'r') as f:
    for line in f:
        if chromRegex.match(line):
            chrom = line.strip()
            chrFaDict[chrom] = list()
        else:
            chrFaDict[chrom].append(line)

os.makedirs(args.output, exist_ok=True)
for chrom in sorted(chrFaDict.keys()):
    if chromKeepRegex.match(chrom):
        realChr = chrom.split('>')[-1]
        chrFaFile = os.path.join(args.output, realChr + '.fa')
        with open(chrFaFile, 'w') as out:
            out.write(chrom + '\n')
            for line in chrFaDict[chrom]:
                out.write(line)
