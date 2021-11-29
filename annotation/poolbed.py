#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
from pybedtools import BedTool
from itertools import combinations

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bed', nargs='+', type=str, required=True,
                    help='input bed6 files')
parser.add_argument('--name', nargs='*', type=str,
                    help='name for input bed (same order with -bed)')
parser.add_argument('--prefix', action='store_true', type=str,
                    default='peak',
                    help='prefix for 4th column of output bed')
parser.add_argument('--strand', action='store_true',
                    default=False,
                    help='overlap bed with strand mode')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='output pooling bed result')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

beds = args.bed

poolBedtoolObj = BedTool()

bedtoolObjList = list(map(lambda x: BedTool(x), beds))

overlapAll = bedtoolObjList[0]
for i in range(1:len(beds)):
    overlapAll = overlapAll.intersect(beds[i], s=args.strand)

