#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################################
#     handle the bed format data         #
#          2018.3.23                     #
##########################################
__author__ = "K.R.Chow"
__version__ = "v1.1"

import sys, os
from functools import reduce
import datetime

############ commands operation ############

# return overlapped-length if bed locus overlapped
def overlap(locusA, locusB):
    # return false if not valid interval
    if (locusA[1] - locusA[0] < 1) or (locusB[1] - locusB[0] < 1):
        return False
    else:
        # overlap length of intervals
        distance = min(locusA[1], locusB[1]) - max(locusA[0], locusB[0])
        return max(0, distance)

# return merged locus
def merge(locusA, locusB, distance=0):
    # distance: Maximum distance between features allowed for intersecting
    # distance = 0 means merge bed like [120,129] and [129,130] to [120, 130]
    overlapLength = overlap(locusA, locusB)
    if overlapLength:
        if overlapLength >= distance:
            return [min(locusA[0], locusB[0]), max(locusA[1], locusB[1])]
        else:
            return False
    else:
        tempDis = min(locusA[1], locusB[1]) - max(locusA[0], locusB[0])
        if tempDis >= distance:
            return [min(locusA[0], locusB[0]), max(locusA[1], locusB[1])]
        else:
            return False

# return intersected locus
def intersect(locusA, locusB, fracA=0.0, fracB=0.0):
    # fracA : overlap length / length of locusA
    # fracB : overlap length / length of locusB
    overlapLength = overlap(locusA, locusB)
    if overlapLength:
        lengthA = locusA[1] - locusA[0]
        lengthB = locusB[1] - locusB[0]
        fracAPercent = overlapLength / lengthA
        fracBPercent = overlapLength / lengthB
        if fracAPercent > float(fracA) and fracBPercent > float(fracB):
            return [max(locusA[0], locusB[0]), min(locusA[1], locusB[1])]
        else:
            return False
    else:
        return False

# format bed12 to exon, intron
def decodeBed12(row):
    start = int(row[1])
    end =  int(row[2])
    thickStart = int(row[6])
    thickEnd = int(row[7])
    blockSizeList = [int(i) for i in row[10].split(',') if i]
    blockStartList = [int(i) for i in row[11].split(',') if i]
    blockList = list(map(lambda x,y:[x + start, x + start + y],
        blockStartList, blockSizeList))
    intronList = list()
    for i in range(len(blockList) - 1):
        intronList.append([blockList[i][1], blockList[i+1][0]])
    if thickStart == thickEnd:
        decodeList = [blockList, intronList]
        return decodeList
    else:
        # decodeList: [[exonblock], [intronblock], [thickup, thick, thickdown]],
        # thickup and thickdown are relative to genome, not strand
        decodeList = [blockList, intronList, [[], [], []]]
        thickStartLocus = [thickStart, thickStart + 1]
        thickEndLocus = [thickEnd - 1, thickEnd]
        thickStartBoolIndex = list(map(lambda x:overlap(thickStartLocus, x),
            blockList)).index(1)
        thickEndBoolIndex = list(map(lambda x:overlap(thickEndLocus, x),
            blockList)).index(1)

        for i in range(len(blockList)):
            blockStart = blockList[i][0]
            blockEnd = blockList[i][1]
            if i < thickStartBoolIndex:
                decodeList[-1][0].append([blockStart, blockEnd])
            elif i == thickStartBoolIndex:
                if thickStart > blockStart:
                    decodeList[-1][0].append([blockStart, thickStart])
                if i == thickEndBoolIndex:
                    decodeList[-1][1].append([thickStart, thickEnd])
                    if thickEnd < blockEnd:
                        decodeList[-1][2].append([thickEnd, blockEnd])
                else:
                    decodeList[-1][1].append([thickStart, blockEnd])
            elif i == thickEndBoolIndex:
                decodeList[-1][1].append([blockStart, thickEnd])
                if thickEnd < blockEnd:
                    decodeList[-1][2].append([thickEnd, blockEnd])
            else:
                decodeList[-1][2].append([blockStart, blockEnd])
    return decodeList

if __name__ == '__main__':
    row = ['chr1','8423769','8424898','ENST00000464367','1000','-','8423770',
        '8423771','0','2','546,93,','0,1036,']
    print(overlap([120,129], [129,130]))
    print(intersect([120,129], [129,130]))
    print(merge([120,129], [129,130], distance=0))
    print(decodeBed12(row))
