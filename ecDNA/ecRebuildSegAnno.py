#!/usr/bin/env python3
import os
import sys
import argparse
import tempfile
import subprocess
import re
from collections import defaultdict
from copy import copy
import operator
from multiprocessing import Pool, Manager
from datetime import datetime

# in-house module
import bedutils

parser = argparse.ArgumentParser(
    description="This script is used for rebuilding segment annotation from annoAaSegBed.py",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-m', '--main', action='store', type=str,
                    required=True,
                    help='annotation output from ecAnnoAaSegBed.py --bed12')
parser.add_argument('-e', '--extra', action='store', type=str,
                    help='annotation output from ecAnnoAaSegBed.py --bed6')
parser.add_argument('--frac', action='store', type=float,
                    default=1,
                    help='cutoff of fraction of genes covered by segment')
parser.add_argument('--cn', action='store', type=float,
                    default=0.05,
                    help='cutoff of copy number')
parser.add_argument('-o', '--output', action='store', type=str,
                    required=True,
                    help='output directory')
parser.add_argument('--otype', action='store', type=str,
                    choices=['whole', 'right', 'left', 'overlay'],
                    help='overlap type of annotation to ecBed')
parser.add_argument('-p', '--prefix', action='store', type=str,
                    required=True,
                    help='prefix of output')
parser.add_argument('--thread', action='store', type=int,
                    default=1,
                    help='thread for running')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def unifyAnno(rebuildExtraDict, discardDict, uniqSegId, chrom, mainList, extraList, repeatPriList):
    discardList = list()
    ## rebuild extra annotation
    reBuildList = list()
    mainLocusList = list(map(lambda x: [chrom, x[8], x[9]], mainList))
    extraList = extraList
    extraList.sort(key = operator.itemgetter(8, 9))
    for row in extraList:
        geneClass = row[4]
        ## element coordinate
        eleLocus = [chrom, row[8], row[9]]
        mainSkipCount = 0
        for i in range(len(mainList)):
            bedops = bedutils.bedops(eleLocus, mainLocusList[i], s=False).intersect()
            if bedops.ibool is True:
                mainGeneClass = mainList[i][4]
                if geneClass in subGeneClassList:
                    reBuildList.append(copy(row))
                else:
                    mainRow = mainList[i]
                    discardRow = copy(row)
                    discardRow.extend(mainRow[0:11])
                    discardList.append(copy(discardRow))
                    break
            else:
                mainSkipCount += 1
        if mainSkipCount == len(mainList):
            if bool(extraList) is False:
                reBuildList.append(copy(row))
            else:
                extraList = copy(reBuildList)
                extraLocusList = list(map(lambda x: [chrom, x[8], x[9]], extraList))
                ## record keep flag number
                appendFlagNum = 0
                for i in range(len(extraList)):
                    preLocus = extraLocusList[i]
                    bedops = bedutils.bedops(eleLocus, preLocus, s=False).intersect()
                    if bedops.ibool is True:
                        preRow = extraList[i]
                        if row[17] == 'repeat' and preRow[17] == 'repeat':
                            typeIndex = repeatPriList.index(row[18])
                            preTypeIndex = repeatPriList.index(preRow[18])
                            if typeIndex < preTypeIndex:
                                ## keep one with higer priority
                                discardRow = copy(preRow)
                                discardRow.extend(row[0:11])
                                discardList.append(copy(discardRow))
                                reBuildList[i] = copy(row)
                                break
                            elif typeIndex == preTypeIndex:
                                ## keep one with more overlapped length
                                preOverLen = int(preRow[30])
                                overLen = int(row[30])
                                if overLen > preOverLen:
                                    discardRow = copy(preRow)
                                    discardRow.extend(row[0:11])
                                    discardList.append(copy(discardRow))
                                    reBuildList[i] = copy(row)
                                break
                            else:
                                discardRow = copy(row)
                                discardRow.extend(preRow[0:11])
                                discardList.append(copy(discardRow))
                                break
                        else:
                            appendFlagNum += 1
                    else:
                        appendFlagNum += 1
                if appendFlagNum == len(extraLocusList):
                    reBuildList.append(copy(row))
    discardDict[uniqSegId] = discardList
    rebuildExtraDict[uniqSegId] = reBuildList

# annotation format
# 0-13: ['#chr', 'start', 'end', 'id', 'copyCount', 'strand', 'sample', 'cycleId', 'cycleLength', 'segUniqApId', 'segLength', 'segNum', 'cycleClass', 'customCycleClass']
# 14-21: ["geneId", "geneName", "synonyms", "description", "geneClass", "geneType", "txId", "txType"]
# 22-33: ["annoStart", "annoEnd", "annoStrand", "ctype", "otype", "orient", "cloverh", "croverh", "absoverh", "overlapLength", "fracA", "fracB"]

startTime = datetime.now()
# built up gene-type pairwise relationships
subGeneClassList = ['intergenic', 'miRNA', 'snoRNA', 'rRNA']
mainAnnoDict = defaultdict(list)
ecBedInforDict = defaultdict(list)
filterUniqSegIdList = list()
headerRow  = list()
discardHeaderRow = list()
with open(args.main, 'r') as f:
    headerRow = f.readline().strip().split('\t')
    discardHeaderRow = copy(headerRow)
    discardHeaderRow.extend(copy(headerRow[0:11]))
    for line in f:
        row = line.strip('\n').split('\t')
        uniqSegId = row[3]
        copyNum = float(row[4])
        otype = row[26]
        fracB = float(row[33])
        geneClass = row[18]
        ## record ecDNA bed
        if uniqSegId not in ecBedInforDict:
            ecBedInforDict[uniqSegId] = row[0:14]
        ## filter by copyNum and overlap fraction
        if copyNum < args.cn:
            continue
        if fracB < args.frac:
            if geneClass != 'intergenic':
                continue
        if bool(args.otype) is True:
            if otype != args.otype:
                continue
        if geneClass in subGeneClassList:
            continue
        ## gene coordinate
        row[22] = int(row[22])
        row[23] = int(row[23])
        mainAnnoDict[uniqSegId].append(copy(row[14:]))
        ## record uniqSegId
        if uniqSegId not in filterUniqSegIdList:
            filterUniqSegIdList.append(uniqSegId)

extraAnnoDict = defaultdict(list)
with open(args.extra, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip('\n').split('\t')
        uniqSegId = row[3]
        copyNum = float(row[4])
        fracB = float(row[33])
        ## record ecDNA bed
        if uniqSegId not in ecBedInforDict:
            ecBedInforDict[uniqSegId] = row[0:14]
        ## filter by copyNum and overlap fraction
        if copyNum < args.cn:
            continue
        if fracB < args.frac:
            if geneClass != 'intergenic':
                continue
        ## gene coordinate
        row[22] = int(row[22])
        row[23] = int(row[23])
        extraAnnoDict[uniqSegId].append(copy(row[14:]))
        ## record uniqSegId
        if uniqSegId not in filterUniqSegIdList:
            filterUniqSegIdList.append(uniqSegId)

discardDict = defaultdict(list)
repeatPriList = ["snRNA", "srpRNA", "scRNA", "LINE", "LTR", "SINE", "RC", "SINE?", "Satellite"] 
repeatPriList += ["Low_complexity", "Simple_repeat", "DNA", "DNA?", "Other", "Unknown"]

if bool(args.thread) is False:
    try:
        thread = int(os.environ['SLURM_JOB_CPUS_PER_NODE'].split('(')[0]) * int(os.environ['SLURM_JOB_NUM_NODES'])
    except KeyError as e:
        thread = os.cpu_count()
    else:
        thread = os.cpu_count()
else:
    thread = args.thread
## run in parallel
pool = Pool(processes=thread)
rebuildExtraDict = Manager().dict()
discardDict = Manager().dict()
filterUniqSegIdList.sort()
for uniqSegId in filterUniqSegIdList:
    chrom = ecBedInforDict[uniqSegId][0]
    ## main annotated genes, ENSG
    mainList = mainAnnoDict[uniqSegId]
    ## repeat, snoRNA...
    extraList = extraAnnoDict[uniqSegId]
    pool.apply_async(unifyAnno, args=(rebuildExtraDict, discardDict, uniqSegId, chrom, mainList, extraList, repeatPriList))
pool.close()
pool.join()

keepFile = os.path.join(args.output, args.prefix + '.rebuild.txt')
with open(keepFile, 'w') as out:
    out.write('\t'.join(headerRow) + '\n')
    allUniqSegIdList = sorted(ecBedInforDict.keys())
    for uniqSegId in allUniqSegIdList:
        bedRow = ecBedInforDict[uniqSegId]
        if uniqSegId in filterUniqSegIdList:
            mainAnnoList = mainAnnoDict[uniqSegId]
            extraAnnoList = rebuildExtraDict[uniqSegId]
            annoList = mainAnnoList + extraAnnoList
            ## sort annotations by coordinates
            annoList.sort(key = operator.itemgetter(8, 9))
            for annoRow in annoList:
                annoRow[8] = str(annoRow[8])
                annoRow[9] = str(annoRow[9])
                row = bedRow + annoRow
                out.write('\t'.join(row) + '\n')
        else:
            annoRow = ['intergenic', 'intergenic', 'na', 'na', 'intergenic', 'intergenic', 'na']
            annoRow += ['na', '.', '.', '.', 'na', 'na', 'na', '0', '0', '0', '0', '0', '0']
            row = bedRow + annoRow
            out.write('\t'.join(row) + '\n')

discardFile = os.path.join(args.output, args.prefix + '.discard.txt')
with open(discardFile, 'w') as out:
    out.write('\t'.join(discardHeaderRow) + '\n')
    for uniqSegId in filterUniqSegIdList:
        if uniqSegId not in discardDict:
            continue
        bedRow = ecBedInforDict[uniqSegId]
        for discardRow in discardDict[uniqSegId]:
            row = bedRow + discardRow
            out.write('\t'.join(map(str,row)) + '\n')

endTime = datetime.now()
print('Duration: {}'.format(endTime - startTime))
