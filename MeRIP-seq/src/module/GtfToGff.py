#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################################
#     designed for formatting            #
#     gff3 and gtf files to bed format   #
#            2017.12.1                   #
##########################################
__author__ = "K.R.Chow"
__version__ = "v1.0"

import sys
import os
import re
from multiBioPro import MultiSys
from collections import defaultdict


#############################################################################
# argument examples:
#   --gff3, ENSEMBL source: ToBed12(file, 'gff3', 'gene', 'mRNA,ncRNA,lnc_RNA,
#                                   pre_miRNA,RNase_MRP_RNA,rRNA,snoRNA,snRNA,SRP_RNA,tRNA', 'ID,Parent', 'list'):
#   --gff3, GENCODE source: ToBed12(file, 'gff3', 'gene', 'transcript', 'gene_id,transcript_id', 'list'):
#############################################################################


saveAccepts = ['list', 'string']


def ToGff3(Object):
    return MainParser(Object)

def MainParser(Object):
    lineList = list()
    #
    if isinstance(Object, str):
        with open(Object, 'r') as f:
            lineList = f.readlines()
    elif isinstance(Object, list):
        lineList = copy.copy(Object)
    else:
        MultiSys.Error(['Unkown object type!'])
    #
    for i in range(len(lineList)):
        if re.search('(^#.*)|(^$)', lineList[i]) is not None:
            lineList[i] = lineList[i].rstrip()
            continue
        lineContList = lineList[i].rstrip().split('\t')
        attributeList = lineContList[-1].replace('"', '').replace(';', '').split(' ')
        tempDict = dict()
        tempDict = MultiSys.List2Dict(attributeList)
        tempList = list()
        #
        if lineContList[2] == 'gene':
            IDVal = tempDict['gene_id']
            del tempDict['gene_id']
            tempList.append('ID=' + IDVal)
        elif lineContList[2] =='transcript':
            IDVal = tempDict['transcript_id']
            parentID = tempDict['gene_id']
            del tempDict['gene_id']
            del tempDict['transcript_id']
            for x in sorted(tempDict.keys()):
                if re.match(r'^gene', x):
                    del tempDict[x]
            tempList.append('ID=' + IDVal)
            tempList.append('Parent=' + parentID)
        else:
            parentID = tempDict['transcript_id']
            del tempDict['gene_id']
            del tempDict['transcript_id']
            for x in sorted(tempDict.keys()):
                if re.match(r'(^gene)|(^transcript)', x):
                    del tempDict[x]
            tempList.append('Parent=' + parentID)
        tempList.extend([x + '=' + tempDict[x] for x in sorted(tempDict.keys())])
        lineContList[-1] = ';'.join(tempList)
        lineList[i] = '\t'.join(lineContList)
    return lineList


def Test():
    file = sys.argv[1]
    testList = ToGff3(file)
    for x in range(len(testList)):
        print(testList[x])


if __name__ == '__main__':
    Test()

