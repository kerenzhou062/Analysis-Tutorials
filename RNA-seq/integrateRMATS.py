#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from copy import copy
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input rMATS results directory')
parser.add_argument('-pval', action='store', type=float,
                    default=1,
                    help='cutoff of pvalue')
parser.add_argument('-fdr', action='store', type=float,
                    default=0.1,
                    help='cutoff of FDR')
parser.add_argument('-diff', action='store', type=float,
                    default=0.1,
                    help='cutoff of IncLevelDifference')
parser.add_argument('-type', action='store', type=str,
                    choices=['JC', 'JCEC'],
                    default='JCEC',
                    help='input type of rMATS results')
parser.add_argument('-grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('-grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output result matrix')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

##public arguments

if bool(args.grepKept):
    kept = True
    regex = re.compile(r'{0}'.format(args.grepKept))
elif bool(args.grepExpel):
    kept = False
    regex = re.compile(r'{0}'.format(args.grepExpel))
else:
    regex = False
rmatsFileList = sorted(glob(os.path.join(args.input, '*.MATS.'+ args.type + '.txt')))

eventFileDict = defaultdict(dict)
for i in range(len(rmatsFileList)):
    rmatsFile = rmatsFileList[i]
    if bool(regex):
        if kept:
            if bool(regex.search(rmatsFile)) is False:
                continue
        else:
            if bool(regex.search(rmatsFile)) is True:
                continue
    rmatsFileName = os.path.basename(rmatsFile)
    event = rmatsFileName.split('.')[0]
    eventFileDict[event] = rmatsFile

geneEventDict = defaultdict(dict)
for event in eventFileDict.keys():
    with open(eventFileDict[event], 'r') as f:
        __ = f.readline()
        for line in f:
            row = line.strip().split('\t')
            pval = float(row[-5])
            fdr = float(row[-4])
            diff = float(row[-1])
            if pval <= args.pval and fdr <= args.fdr and abs(diff) >= args.diff:
                geneId = row[1].replace('"', '')
                geneName = row[2].replace('"', '')
                geneEventDict[geneId]['name'] = geneName
                if 'event' not in geneEventDict[geneId]:
                    geneEventDict[geneId]['event'] = defaultdict(dict)
                    geneEventDict[geneId]['event'][event] = defaultdict(list)
                else:
                    if event not in geneEventDict[geneId]['event']:
                        geneEventDict[geneId]['event'][event] = defaultdict(list)
                if diff < 0:
                    geneEventDict[geneId]['event'][event]['down'].append(diff)
                else:
                    geneEventDict[geneId]['event'][event]['up'].append(diff)

with open(args.output, 'w') as out:
    keyEvent = sorted(eventFileDict.keys())
    eventList = list()
    for event in keyEvent:
        eventList.append(event + '_up')
        eventList.append(event + '_down')
    for event in keyEvent:
        eventList.append(event + '_IncDiff')
    row = ['geneID', 'geneName']
    row.extend(eventList)
    out.write('\t'.join(row) + '\n')
    for geneId in sorted(geneEventDict.keys()):
        geneName = geneEventDict[geneId]['name']
        row = [geneId, geneName]
        for event in keyEvent:
            if event not in geneEventDict[geneId]['event']:
                row.append('0')
                row.append('0')
            else:
                upNum = len(geneEventDict[geneId]['event'][event]['up'])
                downNum = len(geneEventDict[geneId]['event'][event]['down'])
                row.append(str(upNum))
                row.append(str(downNum))
        for event in keyEvent:
            if event not in geneEventDict[geneId]['event']:
                row.append('0')
            else:
                diffList = copy(geneEventDict[geneId]['event'][event]['up'])
                diffList.extend(copy(geneEventDict[geneId]['event'][event]['down']))
                maxDiff = 0
                for diff in diffList:
                    if abs(maxDiff) < abs(diff):
                        maxDiff = diff
                row.append(str(maxDiff))
        out.write('\t'.join(row) + '\n')

