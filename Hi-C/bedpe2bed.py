#!/usr/bin/env python3
import os
import sys
import argparse
import re
import operator
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', action='store', type=str, required=True,
                    help='input bedpe')
parser.add_argument('-p', '--prefix', action='store', type=str,
                    default="bedpe",
                    help='prefix for name column')
parser.add_argument('-s', '--score', action='store', type=int,
                    default=11,
                    help='0-based index column used for score column of bed')
parser.add_argument('-o', '--output', action='store', type=str, required=True,
                    help='the output bed')
parser.add_argument('--addchr', action='store_true',
                    default=False,
                    help='add "chr" to the chromosome name')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

if bool(args.addchr) is True:
    prechr = "chr"
else:
    prechr = ""

regex = re.compile(r'^#')
lcount = 0
bedrowList = list()
with open(args.input, 'r') as f:
    for line in f:
        lcount += 1
        name = args.prefix + "|" + str(lcount)
        if regex.match(line):
            continue
        else:
            row = line.strip().split('\t')
            chr1 = prechr + row[0]
            start1 = row[1]
            end1 = row[2]
            name1 = name + "|1"
            strand1 = row[8]
            chr2 = prechr + row[3]
            start2 = row[4]
            end2 = row[5]
            name2 = name + "|2"
            strand2 = row[9]
            score = row[args.score]
            bedrowList.append([chr1, start1, end1, name1, score, strand1])
            bedrowList.append([chr2, start2, end2, name2, score, strand2])

bedrowList.sort(key = operator.itemgetter(0, 1))
with open(args.output, 'w') as output:
    for bedrow in bedrowList:
        output.write('\t'.join(map(str, bedrow)) + '\n')
