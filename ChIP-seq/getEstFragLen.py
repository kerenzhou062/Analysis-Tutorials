#!/usr/bin/env python3
import os
import sys
import argparse
import re
import pathlib

#usage: runChipPeakCallBash.py or runChipPeakCallBash.py <fastq dir>

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-extsize', action='store', type=int,
                    help='--extsize parameter (macs2, histone:147)')
parser.add_argument('-input', action='store', type=str,
                    default='./', help='the root IP/input directory (eg. alignment/HepG2_CTCF_WT/IP)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

# public arguments
extsize = args.extsize
basepath = os.path.realpath(args.input)

if bool(args.extsize) is True:
    extsize = args.extsize
else:
    qcFileList = pathlib.Path(basepath).rglob("*.trim.filt.sample.15.tagAlign.gz.cc.qc")
    valideFileNum = 0
    estFragLenSum = 0
    for qcFile in qcFileList:
        qcFileName = qcFile.name
        if re.search(r'_IP_', qcFileName):
            valideFileNum += 1
            with qcFile.open() as f:
                line = f.readline()
                row = line.strip().split('\t')
                estFragLen = int(row[2])
                estFragLenSum += estFragLen
    if estFragLenSum == 0:
        extsize = 147
    else:
        extsize = int(estFragLenSum / valideFileNum)
print(extsize)
