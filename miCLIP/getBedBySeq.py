#!/usr/bin/env python3
import os
import sys
import argparse
import tempfile
from collections import defaultdict
import subprocess

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-fasta', action='store', type=str, required=True,
                    help='genome fasta')
parser.add_argument('-bed', action='store', type=str, required=True,
                    help='input bed file (bed6+)')
parser.add_argument('-seq', action='store', type=str, required=True,
                    help='seqeuence to be kept (case insensitive)')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output bed file')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

bed6Tmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
skipLineRow = []
keyBedDict = defaultdict(list)
with open(args.bed, 'r') as f, open(bed6Tmp.name, 'w') as temp:
    for line in f:
        row = line.strip().split('\t')
        if line[0] == '#':
            skipLineRow.append(line)
        else:
            chrom = row[0]
            start = row[1]
            end = row[2]
            strand = row[5]
            key = "{0}:{1}-{2}({3})".format(chrom, start, end, strand)
            keyBedDict[key].append(line)
            bed6Line = '\t'.join(row[:6]) + '\n'
            temp.write(bed6Line)

seq = args.seq.upper()

getFastaTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
command = "bedtools getfasta -fi {0} -bed {1} -tab -s > {2}".format(args.fasta, bed6Tmp.name, getFastaTmp.name)
subprocess.call(command, shell=True)
keptKeyList = []
with open(getFastaTmp.name, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        key = row[0]
        keySeq = row[1].upper()
        if seq == keySeq:
            if key not in keptKeyList:
                keptKeyList.append(key)

with open(args.output, 'w') as out:
    if len(skipLineRow) != 0:
        for skipLine in skipLineRow:
            out.write(skipLine)
    for key in keptKeyList:
        for line in keyBedDict[key]:
            out.write(line)
