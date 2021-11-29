#!/usr/bin/env python3
import os
import sys
import argparse
from glob import glob
from copy import copy
import re
import subprocess
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="This script is used for converging ecDNA results from runPipeAA.sh to bed6+",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', action='store', type=str, required=True,
                    help='directory for storing ecDNA results')
parser.add_argument('--genome', action='store', type=str, choices=['hg38', 'hg19'],
                    default='hg19',
                    help='input genome build')
parser.add_argument('--grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('--grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('--length', action='store', type=int,
                    default=10000,
                    help='length in bp to be recognized as big circular fragment')
parser.add_argument('-o', '--output', action='store', type=str, required=True,
                    help='output bed')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def CycleClassifyCustom(vertex, segDict, trueSegIndexList, segIdList, segStrandList):
    # classify cicyles
    # ref. https://www.nature.com/articles/s41586-019-1913-9
    allSegNum = len(trueSegIndexList)
    cycleClass = 'other'
    chromDict = defaultdict(int)
    segSetDict = defaultdict(dict)
    trueSegIdList = list()
    trueSegStrandList = list()
    for i in trueSegIndexList:
        segId = segIdList[i]
        strand = segStrandList[i]
        trueSegIdList.append(segId)
        trueSegStrandList.append(strand)
        segChr, segStart, segEnd = segDict[segId]
        chromDict[segChr] += 1
        if segId not in segSetDict:
            segSetDict[segId] = defaultdict(list)
        segSetDict[segId]['strand'].append(strand)
        segSetDict[segId]['order'].append(i)
    chromNum = len(chromDict.keys())
    segUniqNum = len(segSetDict.keys())
    ## decode seg composition
    segPosDict = defaultdict(int)
    strandList = list()
    for segId in sorted(segSetDict.keys()):
        eSegStrandList = segSetDict[segId]['strand']
        eSegOrderList = segSetDict[segId]['order']
        strandList.extend(eSegStrandList)
        eSegstrandNum = len(set(eSegStrandList))
        eSegOrderNum = len(eSegOrderList)
        linkSeg = 0
        skipSeg = 0
        if eSegOrderNum > 1:
            for i in range(1,eSegOrderNum):
                gap = eSegOrderList[i] - eSegOrderList[i-1]
                if gap == 1:
                    linkSeg += 1
                else:
                    skipSeg += 1
            if linkSeg > 0 and skipSeg == 0:
                if eSegstrandNum == 1:
                    segPosDict['link'] += 1
                else:
                    segPosDict['link_inversion'] += 1
            else:
                if eSegstrandNum == 1:
                    segPosDict['skip'] += 1
                else:
                    segPosDict['skip_inversion'] += 1
    strandNum = len(set(strandList))
    ## to classify
    if vertex == 'circular':
        cycleClass = 'circular'
    elif allSegNum == 1:
        cycleClass = 'single'
    else:
        if chromNum > 2:
            if strandNum > 1:
                cycleClass = 'chromoplexy_inversion'
            else:
                cycleClass = 'chromoplexy'
        else:
            if allSegNum == segUniqNum:
                if chromNum == 1:
                    if strandNum == 1:
                        cycleClass = 'deletion'
                    else:
                        if segUniqNum == 3:
                            if trueSegStrandList[0] == trueSegStrandList[2]:
                                cycleClass = 'reciprocal_inversion'
                            else:
                                cycleClass = 'inversion'
                        else:
                            if segUniqNum < 5:
                                cycleClass = 'inversion'
                            else:
                                cycleClass = 'chromothripsis'
                else:
                    cycleClass = 'translocation'
            else:
                link = segPosDict['link']
                linkInver = segPosDict['link_inversion']
                skip = segPosDict['skip']
                skipInver = segPosDict['skip_inversion']
                if chromNum == 1:
                    if segUniqNum == 1 and allSegNum == 2 and trueSegIdList[0] == trueSegIdList[1]:
                        cycleClass = 'self_inversion'
                    elif link == 1 and linkInver == 0 and skip == 0 and skipInver == 0:
                        cycleClass = 'tandem_duplication'
                    elif link == 0 and linkInver == 0 and skip == 1 and skipInver == 0:
                        cycleClass = 'local_distant_cluster'
                    elif skip > 1:
                        cycleClass = 'local_n_jump'
                    elif link == 0 and linkInver == 0 and skip == 0 and skipInver > 0:
                        cycleClass = 'jump_inversion'
                    elif link == 0 and linkInver > 0 and skip == 0 and skipInver == 0:
                        cycleClass = 'tandem_duplication_inversion'
                    ## foldback_inversions
                    if segUniqNum == 2 and allSegNum == 3:
                        if strandNum > 1:
                            if trueSegIdList[0] == trueSegIdList[2] and trueSegIdList[0] != trueSegIdList[1]:
                                if trueSegStrandList[0] == '+' and trueSegStrandList[2] == '-':
                                    cycleClass = 'foldback_inversion'
                                elif trueSegStrandList[0] == '-'and trueSegStrandList[2] == '+':
                                    cycleClass = 'foldback_inversion'
                else:
                    if link == 0 and linkInver == 0 and skip == 1 and skipInver == 0:
                        cycleClass = 'local_distant_cluster'
                    else:
                        cycleClass = 'complexity'
    return cycleClass

# main program
if bool(args.grepKept):
    kept = True
    regex = re.compile(r'{0}'.format(args.grepKept))
elif bool(args.grepExpel):
    kept = False
    regex = re.compile(r'{0}'.format(args.grepExpel))
else:
    regex = False

if args.genome == 'hg38':
    args.genome = 'GRCh38'

cyclesFileList = sorted(glob(os.path.join(args.input, '**', '*_cycles.txt'), recursive=True))
graphFileList = sorted(glob(os.path.join(args.input, '**', '*_graph.txt'), recursive=True))

cyclesFiles = list()
graphFileDict = defaultdict(dict)
for cyclesFile in cyclesFileList:
    if bool(regex):
        if kept:
            if bool(regex.search(cyclesFile)) is False:
                continue
        else:
            if bool(regex.search(cyclesFile)) is True:
                continue
    cyclesFiles.append(cyclesFile)

for graphFile in graphFileList:
    if bool(regex):
        if kept:
            if bool(regex.search(graphFile)) is False:
                continue
        else:
            if bool(regex.search(graphFile)) is True:
                continue
    graphFileName = os.path.basename(graphFile)
    graphFileDict[graphFileName] = graphFile

if len(cyclesFiles) != len(graphFileDict.keys()):
    print('Number of cycles files is not equal to graph files')
    sys.exit()

for cyclesFile in cyclesFiles:
    # eg, K562_encode_amplicon55_cycles.txt
    # eg, K562_encode_amplicon55_graph.txt
    dirName = os.path.dirname(cyclesFile)
    cyclesFileName = os.path.basename(cyclesFile)
    nameRow = cyclesFileName.split('_')
    nameLen = len(nameRow)
    graphFileName = '_'.join(nameRow[0:nameLen-1] + ['graph.txt'])
    if graphFileName not in graphFileDict:
        print(cyclesFileName + 'did not have corresponding graph file')
        sys.exit()

# convergeDict: 'sample'->'ampliconId'->cycleDict
convergeDict = defaultdict(dict)
cycleLineDict = defaultdict(dict)
for cyclesFile in cyclesFiles:
    # eg, K562_encode_amplicon55_cycles.txt
    # eg, K562_encode_amplicon55_graph.txt
    dirName = os.path.dirname(cyclesFile)
    cyclesFileName = os.path.basename(cyclesFile)
    nameRow = cyclesFileName.split('_')
    nameLen = len(nameRow)
    graphFileName = '_'.join(nameRow[0:nameLen-1] + ['graph.txt'])
    # get graph file
    if graphFileName not in graphFileDict:
        print(cyclesFileName + 'did not have corresponding graph file')
        sys.exit()
    graphFile = graphFileDict[graphFileName]
    # get amplicon information
    ampliconId = nameRow[nameLen - 2]
    sample = '_'.join(nameRow[0:nameLen-2])
    uniqAmpliconId = '_'.join([sample, ampliconId])
    ## cycleDict:'seg'->[pos]; 'cycle'->'id'->'seg','copyCount', 'seg' = [segIdList, segStrandList]
    cycleDict = defaultdict(dict)
    ## deal with segment
    breakpointDict = defaultdict(dict)
    with open(cyclesFile, 'r') as f:
        for line in f:
            row = re.split('\t|;|=', line.strip('\n'))
            if row[0] == 'Segment':
                segId = row[1]
                segChr = row[2]
                # 0-based
                segStart = int(row[3])
                segEnd = int(row[4])
                if segStart == segEnd:
                    breakpointDict[segId] = True
                else:
                    cycleDict['seg'][segId] = [segChr, segStart, segEnd]
    ## deal with cycles
    with open(cyclesFile, 'r') as f:
        for line in f:
            row = re.split('\t|;|=', line.strip('\n'))
            if row[0] == 'Cycle':
                # record cycle line and remove duplicate rerords
                uniqCycleLine = '_'.join([uniqAmpliconId, re.sub(r'Cycle=\d+;', '', line.strip('\n'))])
                if uniqCycleLine not in cycleLineDict:
                    cycleLineDict[uniqCycleLine] = 1
                else:
                    continue
                cycleId = row[1]
                copyCount = row[3]
                cycleSegRow = row[-1].split(',')
                segIdList = []
                segStrandList = []
                for cycleSeg in cycleSegRow:
                    cycleSeg = cycleSeg.replace('+', ':+')
                    cycleSeg = cycleSeg.replace('-', ':-')
                    segRow = cycleSeg.split(':')
                    if len(segRow) != 2:
                        continue
                    else:
                        if segRow[0] not in breakpointDict:
                            segIdList.append(segRow[0])
                            segStrandList.append(segRow[1])
                    if cycleId not in cycleDict['cycle']:
                        cycleDict['cycle'][cycleId] = defaultdict(dict)
                    cycleDict['cycle'][cycleId]['seg'] = [segIdList, segStrandList]
                    cycleDict['cycle'][cycleId]['copyCount'] = copyCount
            else:
                continue
    ## deal with cycle class with AmpliconClassifier
    ampliconAnnoCommand = 'amplicon_annotate_cycle.py -c {0} -g {1} --ref {2} -o ./ --stdout'.format(cyclesFile, graphFile, args.genome)
    ampliconAnnoRes = bytes.decode(subprocess.check_output(ampliconAnnoCommand, shell=True)).split('\n')
    for line in ampliconAnnoRes:
        if bool(line) is False:
            continue
        row = line.strip().split('\t')
        if bool(re.match(r'^Cycle=.+IsCyclicPath=', row[0])) is True:
            annoRow = row[0].split(';')
            cycleId = annoRow[0].split('=')[-1]
            copyCount = float(annoRow[1].split('=')[-1])
            cycleLength = float(annoRow[2].split('=')[-1])
            cyclic = annoRow[3].split('=')[-1]
            cycleClass = annoRow[4].split('=')[-1]
            if cycleClass == "Invalid" and cyclic == "True":
                if cycleLength < args.length:
                    cycleClass = "Invalid"
                else:
                    cycleClass = "Circular"
            else:
                if cycleClass == "Invalid" and cycleLength >= 500:
                    cycleClass = "No_fSCNA"
                else:
                    cycleClass == "Invalid"
            cycleDict['cycle'][cycleId]['CycleClass'] = cycleClass
    convergeDict[sample][ampliconId] = cycleDict

with open(args.output, 'w') as out:
    row = ['#chr', 'start', 'end', 'id', 'copyCount', 'strand', 'sample', \
        'cycleId', 'cycleLength', 'segUniqApId', 'segLength', 'segNum', 'cycleClass', 'customCycleClass']
    out.write('\t'.join(row) + '\n')
    for sample in sorted(convergeDict.keys()):
        for ampliconId in sorted(convergeDict[sample].keys()):
            cycleDict = convergeDict[sample][ampliconId]
            for cycleId in sorted(cycleDict['cycle'].keys()):
                cycleClass = cycleDict['cycle'][cycleId]['CycleClass']
                if bool(cycleClass) is False or cycleClass == 'Invalid':
                    continue
                segIdList, segStrandList = cycleDict['cycle'][cycleId]['seg']
                copyCount = cycleDict['cycle'][cycleId]['copyCount']
                cycleLength = cycleDict['cycle'][cycleId]['length']
                trueSegIndexList = list()
                vertexList = list()
                segNum = len(segIdList)
                for i in range(segNum):
                    if segIdList[i] =='0' :
                        vertexList.append(segIdList[i] + segStrandList[i])
                    else:
                        trueSegIndexList.append(i)
                vertex = ','.join(vertexList)
                if vertex == '0-,0+':
                    vertex = '0+,0-'
                if vertex == '':
                    vertex = 'circular'
                cycleLength = sum(map(lambda x: cycleDict['seg'][segIdList[x]][2] - cycleDict['seg'][segIdList[x]][1], trueSegIndexList))
                trueSegNum = len(trueSegIndexList)
                cycleUniqId = '|'.join([sample, ampliconId, 'cycle-' + cycleId])
                customCycleClass = CycleClassifyCustom(vertex, cycleDict['seg'], trueSegIndexList, segIdList, segStrandList)
                if customCycleClass == 'single':
                    continue
                for i in range(trueSegNum):
                    segIndex = trueSegIndexList[i]
                    segId = segIdList[segIndex]
                    strand = segStrandList[segIndex]
                    segChr, segStart, segEnd = cycleDict['seg'][segId]
                    uniqSegId = '|'.join([sample, ampliconId, 'cycle-' + cycleId, 'seg-' + str(i+1)])
                    segUniqApId = 'seg-' + segId
                    segLength = segEnd - segStart
                    row = [segChr, str(segStart), str(segEnd), uniqSegId, copyCount, strand, sample, \
                        cycleUniqId, str(cycleLength), segUniqApId, str(segLength), str(trueSegNum), cycleClass, customCycleClass]
                    out.write('\t'.join(row) + '\n')
