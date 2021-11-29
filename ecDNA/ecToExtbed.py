#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
from copy import copy
import re
# in-house module
import bedutils

parser = argparse.ArgumentParser(
    description="This script is used for transforming coordinate from \
    .rebuild.txt to circle-based bed and fasta",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-e', '--extend', action='store', type=int,
                    default=2000,
                    help='increase the extended size around junction sites')
parser.add_argument('-i', '--input', action='store', type=str,
                    required=True,
                    help='bed output from ecConvergeAaToSegBed.py')
parser.add_argument('-l', '--length', action='store', type=int,
                    default=1000,
                    help='cutoff of cycle length')
parser.add_argument('-s', '--size', action='store', type=str,
                    required=True,
                    help='genome size file')
parser.add_argument('-o', '--output', action='store', type=str,
                    required=True,
                    help='output folder')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def chunkString(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

chrSizeDict = defaultdict(int)
with open(args.size, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        chrom = row[0]
        chromSize = int(row[1])
        chrSizeDict[chrom] = chromSize

sampleCoorDict = defaultdict(dict)
with open(args.input, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        chrom = row[0]
        start = row[1]
        end = row[2]
        segId = row[3]
        score = row[4]
        strand = row[5]
        sample = row[6]
        cycleType = row[12]
        cycleTypeCustom = row[13]
        cycleLength = int(row[8])
        if cycleLength < args.length:
            continue
        coordinate = ';'.join([chrom, start, end, strand])
        if coordinate not in sampleCoorDict[sample]:
            sampleCoorDict[sample][coordinate] = defaultdict(list)
        sampleCoorDict[sample][coordinate]['segId'].append(segId)
        sampleCoorDict[sample][coordinate]['cycleType'].append(cycleType)
        sampleCoorDict[sample][coordinate]['cycleTypeCustom'].append(cycleTypeCustom)

os.makedirs(args.output, exist_ok=True)
for sample in sorted(sampleCoorDict.keys()):
    extendBed = os.path.join(args.output, sample + '.extend.bed')
    coorList = sorted(sampleCoorDict[sample].keys())
    with open(extendBed, 'w') as out:
        for i in range(len(coorList)):
            coordinate = coorList[i]
            chrom, start, end, strand = coordinate.split(';')
            score = str(len(sampleCoorDict[sample][coordinate]['segId']))
            name = 'segment-' + str(i+1)
            chromSize = chrSizeDict[chrom]
            extStart = int(start)- args.extend
            extEnd = int(end) + args.extend
            if extStart < 0:
                extStart = 0
            if extEnd > chromSize:
                extEnd = chromSize
            row = [chrom, str(extStart), str(extEnd), name, score, strand]
            segIds = ','.join(sampleCoorDict[sample][coordinate]['segId'])
            cycleTypes = ','.join(sampleCoorDict[sample][coordinate]['cycleType'])
            cycleTypesCustom = ','.join(sampleCoorDict[sample][coordinate]['cycleTypeCustom'])
            row.extend([str(start), str(end), segIds, cycleTypes, cycleTypesCustom])
            out.write('\t'.join(row) + '\n')
