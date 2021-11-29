#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
from glob import glob
import subprocess
import tempfile
## custom modules
import bedutils

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='This script is used for trimming fasta based on avrage length of sequences in --ref fasta.\
     (the length of sequences in --input must be the same and larger than avrage length)')
parser.add_argument('-r', '--ref', action='store', type=str,
                    required=True,
                    help='genome fasta file')
parser.add_argument('-i','--input', action='store', type=str,
                    required=True,
                    help='input bed file in bed6 or bed12 format')
parser.add_argument('-o', '--output', action='store', type=str,
                    required=True,
                    help='The output fasta')
parser.add_argument('-t','--tag', action='store', type=str,
                    default="background",
                    help='the string attached to the fastaid')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

regex = re.compile(r'^>')
refDict = defaultdict(list)
with open(args.ref, 'r') as f:
    for line in f:
        line = line.strip()
        if regex.match(line):
            fastaid = line
        else:
            refDict[fastaid].append(line)

fastaDict = defaultdict(list)
with open(args.input, 'r') as f:
    for line in f:
        line = line.strip()
        if regex.match(line):
            fastaid = line
        else:
            fastaDict[fastaid].append(line)

seqLengthSum = 0
for fastaid in sorted(refDict.keys()):
    seq = ''.join(refDict[fastaid])
    seqLengthSum += len(seq)

refHalfLenAve = int(seqLengthSum / len(sorted(refDict.keys())) / 2)

with open(args.output, 'w') as out:
    for fastaid in sorted(fastaDict.keys()):
        seq = ''.join(fastaDict[fastaid])
        ## trim fasta
        center = int(len(seq) / 2)
        start = center - refHalfLenAve
        end = center + refHalfLenAve
        finalSeq = seq[start:end]
        fastaid = '|'.join([fastaid, args.tag])
        out.write(fastaid + '\n')
        out.write(finalSeq + '\n')
