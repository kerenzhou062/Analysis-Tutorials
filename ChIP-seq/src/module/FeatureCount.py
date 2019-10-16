#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################
#     designed for calculating expression    #
#     level from featureCount results        #
#            2017.12.12                      #
##############################################
__author__ = "K.R.Chow"

import sys
import os
import argparse
import re
from multiBioPro import MultiSys
from collections import defaultdict

#############################################################################
# argument examples:
#   RpmC(file, mode='rpm', cutoff=0, precision=3):
#############################################################################

thousand = 1e3
milion = 1e6

def Exp(file, mode='rpm', cutoff=0, precision=3):
    if mode == 'rpm':
        return RpmC(file, cutoff, precision)
    elif mode == 'fpkm' or mode == 'rpkm':
        return RpkmC(file, cutoff, precision)
    elif mode == 'tpm':
        return TpmC(file, cutoff, precision)


def RpmC(file, cutoff, precision):
    totalList, lineList = MainParser(file, mode='rpm')
    for i in range(1, len(lineList)):
        countsList = lineList[i][6:]
        if sum(countsList) < cutoff:
            continue
        lineList[i][6:] = list(map(lambda x, y: round(
            x / y, precision), countsList, totalList))
    return lineList


def RpkmC(file, cutoff, precision):
    totalList, lineList = MainParser(file, mode='rpkm')
    for i in range(1, len(lineList)):
        countsList = lineList[i][6:]
        length = lineList[i][5] / thousand
        if sum(countsList) < cutoff:
            continue
        lineList[i][6:] = list(map(lambda x, y: round(
            x / (y * length), precision), countsList, totalList))
    return lineList


def TpmC(file, cutoff, precision):
    totalList, lineList = MainParser(file, mode='tpm')
    for i in range(1, len(lineList)):
        countsList = lineList[i][6:]
        if sum(countsList) < cutoff:
            continue
        lineList[i][6:] = list(map(lambda x, y: round(x / y, precision), countsList, totalList))
    return lineList


def MainParser(file, mode):
    lineList = list()
    totalList = list()
    sampleSize = 0
    with open(file, 'r') as f:
        for line in f.readlines():
            if re.match(r'^#', line):
                continue
            row = line.rstrip('\n').split('\t')
            if row[0] == 'Geneid':
                sampleSize = len(row) - 6
                totalList = [[] for _ in range(sampleSize)]
                lineList.append(row)
            else:
                length = row[5] = int(row[5])
                countsList = list(map(int, (row[6:])))
                if mode == 'tpm':
                    countsList = list(map(lambda x: (x * thousand) / length, countsList))
                    for j in range(sampleSize):
                        totalList[j].append(countsList[j])
                lineList.append(row[0:6] + countsList)
    if mode == 'tpm':
        totalList = list(map(lambda x: x/milion, map(sum, totalList)))
    else:
        sumaryFile = file + '.summary'
        MultiSys.FileExist(sumaryFile)
        with open(sumaryFile, 'r') as f:
            __ = f.readline()
            for line in f.readlines():
                for j in range(sampleSize):
                    readList = line.rstrip('\n').split('\t')[1:]
                    totalList[j].append(int(readList[j]))
        totalList = list(map(lambda x:x/milion, map(sum, totalList)))
    return totalList, lineList
