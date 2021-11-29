#!/usr/bin/env python3
import os
import sys
import argparse
from glob import glob
from copy import copy
import re
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="This script is used for doing the statistis of BWA mapping",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', action='store', type=str, required=True,
                    help='directory for storing ecDNA results')
parser.add_argument('-o', '--output', action='store', type=str, required=True,
                    help='output bed')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

bamStatsFiles = sorted(glob(os.path.join(args.input, '**', 'bam.stats.txt'), recursive=True))

statsDict = defaultdict(dict)
for bamStatsFile in bamStatsFiles:
    sampleName = os.path.split(os.path.split(bamStatsFile)[0])[-1]
    with open(bamStatsFile, 'r') as f:
        totalReadNum = f.readline().strip().split(' ')[0]
        