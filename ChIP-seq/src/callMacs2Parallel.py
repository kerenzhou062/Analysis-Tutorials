#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "K.R.Chow"

import sys
import os
import argparse
import datetime
import shutil
from multiprocessing import Pool, Manager
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='cores for multiprocessing')
parser.add_argument('-dtype', action='store', type=str,
                    default='narrowPeak', help='The --input-file-type \
                    (narrowPeak|broadPeak|bed|gff)')
parser.add_argument('-ip', nargs='+', type=str,
                    help='The ip files')
parser.add_argument('-input', nargs='+', type=str,
                    help='The input file strings')
parser.add_argument('-gsize', action='store', type=str,
                    default='hs', help='genome size')
parser.add_argument('-format', action='store', type=str,
                    help='input format (BED, BAM...)')
parser.add_argument('-shift', action='store', type=int,
                    default=0, help='shift value')
parser.add_argument('-extsize', action='store', type=int,
                    default=147, help='extsize value')
parser.add_argument('-name', nargs='+', type=str,
                    help='The prefix name')
parser.add_argument('-output', nargs='+', type=str,
                    help='The output result directories')
parser.add_argument('-pval', action='store', type=float,
                    default=1e-2, help='The p-value cutoff')
parser.add_argument('-other', action='store', type=str,
                    help='The other commands, "," separated (eg.--broad,--broad)')
parser.add_argument('-rank', action='store', type=str,
                    default='signal.value', help='Which column to use to rank peaks \
                    (idr, signal.value|p.value|q.value)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def runMacs2(ip, control, name, other, output, extsize, gsize, shift, pval, dFormat, dtype, rank):
    # calling peaks
    macs2Command = 'macs2 callpeak -t {ip} \
        -c {control} -f {dFormat} -g {gsize} \
        -p {pval} --nomodel --shift {shift} \
        --extsize {extsize} -n {name} --outdir {output} \
        {other} > {output}/../{name}.macs2.log 2>&1'.format(**vars())
    subprocess.run(macs2Command, shell=True)
    # sort peaks by signal.value | -log10(p-value) | -log10(p-value)
    peakFile = os.path.join(output, name + '_peaks.{0}'.format(dtype))
    broadFlag = False
    if 'broadPeak' in dtype:
        broadFlag = True
    if rank == 'signal.value':
        sortCommand = 'sort -k7,7n {0} -o {0}'.format(peakFile)
    elif rank == 'p.value':
        sortCommand = 'sort -k8,8nr {0} -o {0}'.format(peakFile)
    else:
        sortCommand = 'sort -k9,9nr {0} -o {0}'.format(peakFile)
    subprocess.run(sortCommand, shell=True)
    return 1

# public variables
ipList = args.ip
inputList = args.input
nameList = args.name
otherList = args.other.split(',')
outputList = args.output
extsize = args.extsize
gsize = args.gsize
shift = args.shift
pval = args.pval
dFormat = args.format
dtype = args.dtype
rank = args.rank

# judge elements in same length
if len(ipList) == len(inputList) and len(ipList) == len(nameList):
    if len(ipList) == len(otherList):
        if len(ipList) == len(outputList):
            pass
        else:
            sys.exit('not equal parameters')
    else:
        sys.exit('not equal parameters')
else:
    sys.exit('not equal parameters')

starttime = datetime.datetime.now()
sys.stderr.write("Parallel peak-calling with macs2!\n")

# multiple call macs2
pool = Pool(processes=args.cpu)
for i in range(len(ipList)):
    ip = ipList[i]
    control = inputList[i]
    name = nameList[i]
    other = otherList[i]
    output = outputList[i]
    # try to make dirs
    try:
        shutil.rmtree(output)
    except FileNotFoundError as e:
        pass
    os.makedirs(output, exist_ok=True)
    #runMacs2(ip, control, name, other, output, extsize, gsize, shift, pval, dFormat)
    result = pool.apply_async(runMacs2, args=(ip, control, name, other, 
        output, extsize, gsize, shift, pval, dFormat, dtype, rank))
pool.close()
pool.join()

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("All jobs done!")
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
