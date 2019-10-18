#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "K.R.Chow"

import sys
import os
import argparse
import datetime
import shutil
import math
from multiprocessing import Pool
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-blacklist', action='store', type=str,
                    help='The blacklist bed (from ENCODE)')
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='cores for multiprocessing')
parser.add_argument('-dtype', action='store', type=str,
                    default='narrowPeak', help='The --input-file-type \
                    (narrowPeak|broadPeak|bed|gff)')
parser.add_argument('-output', nargs='+', type=str,
                    help='The output directories')
parser.add_argument('-peaklist', nargs='+', type=str,
                    help='The pooled peaks (not necessary)')
parser.add_argument('-prefix', nargs='+', type=str,
                    help='The output prefixes')
parser.add_argument('-rank', action='store', type=str,
                    default='signal.value', help='Which column to use to rank peaks \
                    (signal.value|p.value|q.value|columnIndex)')
parser.add_argument('-sample', action='store', type=str,
                    help='The input replicate peaks string \
                    ("rep1 rep2,pr_rep1 pr_rep2")')
parser.add_argument('-thresh', action='store', type=float,
                    default=0.05, help='The --soft-idr-threshold')


args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def runIdr(sample, peaklist, prefix, output, blacklist, dtype, thresh, rank):
    logThresh = -math.log(thresh, 10)
    thresh = str(thresh)
    #run IDR on samples and peaklist
    peakListParam = '--peak-list {0}'.format(peaklist) if bool(peaklist) else ''
    idrOutput = os.path.join(output, prefix + '.IDR.' + dtype)
    idrCallCommand = 'idr --samples {sample} {peakListParam} \
        --output-file {idrOutput} --rank {rank} \
        --soft-idr-threshold {thresh} --plot \
        --use-best-multisummit-IDR \
        > {output}/../{prefix}.IDR.Call.log 2>&1'.format(**vars())
    subprocess.run(idrCallCommand, shell=True)

    #filter with IDR --soft-idr-thresh
    idrThredOutput = os.path.join(output, prefix + 
        'IDR.' + thresh + '.' + dtype)
    idrThredCommand = "awk 'BEGIN{{OFS=\"\\t\"}} $12>='{logThresh}' \
        {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {idrOutput} \
        | sort | uniq | sort -k7n,7n > {idrThredOutput}".format(**vars())
    subprocess.run(idrThredCommand, shell=True)

    #filter with blacklist bed
    idrFiltOutput = os.path.join(output, prefix + 
        'IDR.' + thresh + '.filt.' + dtype)
    idrFiltCommand = "bedtools intersect -v -a {idrThredOutput} -b {blacklist} \
        | grep -P 'chr[\\dXY]+[ \t]' | awk 'BEGIN{{OFS=\"\\t\"}} \
        {{if ($5>1000) $5=1000; print $0}}' \
        > {idrFiltOutput}".format(**vars())
    subprocess.run(idrFiltCommand, shell=True)
    return 1

# public variables
sampleList = args.sample.split(',')
peaklistList = args.peaklist
prefixList = args.prefix
outputList = args.output
blacklist = args.blacklist
dtype = args.dtype
thresh = args.thresh
rank = args.rank

# judge elements in same length
if len(sampleList) == len(prefixList) and len(sampleList) == len(outputList):
    if bool(peaklistList):
        if len(sampleList) == len(peaklistList):
            pass
        else:
            sys.exit('not equal parameters')
    else:
        pass
else:
    sys.exit('not equal parameters')

starttime = datetime.datetime.now()
sys.stderr.write("Parallel IDR-running with idr!\n")

# multiple call idr
pool = Pool(processes=args.cpu)
for i in range(len(sampleList)):
    sample = sampleList[i]
    if bool(peaklistList):
        peaklist = peaklistList[i]
    else:
        peaklist = ""
    prefix = prefixList[i]
    output = outputList[i]
    # try to make dirs
    try:
        shutil.rmtree(output)
    except FileNotFoundError as e:
        pass
    os.makedirs(output, exist_ok=True)
    #runIdr(sample, peaklist, prefix, output, blacklist, dtype, thresh, rank)
    result = pool.apply_async(runIdr, args=(sample, peaklist, prefix, output, 
        blacklist, dtype, thresh, rank))
pool.close()
pool.join()

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("All jobs done!")
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
