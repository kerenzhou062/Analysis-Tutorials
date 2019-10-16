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
import copy
from multiBioPro import MultiSys
from collections import defaultdict


#############################################################################
# argument examples:
#   --gff3, ENSEMBL source: ToBed12(file, 'gff3', 'gene', 'mRNA,ncRNA,lnc_RNA,
#                                   pre_miRNA,RNase_MRP_RNA,rRNA,snoRNA,snRNA,SRP_RNA,tRNA', 'ID,Parent', 'list'):
#   --gff3, GENCODE source: ToBed12(file, 'gff3', 'gene', 'transcript', 'gene_id,transcript_id', 'list'):
#############################################################################


saveAccepts = ['list', 'string']


def ToBed12(Object, formatType='gff3', genePattern='gene', txPattern='transcript', identifier='ID,Parent', save='list'):
    MultiSys.AcceptArgs(save, *saveAccepts)
    geneDict, txDict = MainParser(
        Object, formatType, genePattern, txPattern, identifier)
    return DictParser(mode='bed12', feature='tx', save=save, **txDict)


def ToBed6(Object, formatType='gff3', feature='gene', genePattern='gene', txPattern='transcript', identifier='ID,Parent', save='list'):
    MultiSys.AcceptArgs(save, *saveAccepts)
    geneDict, txDict = MainParser(
        Object, formatType, genePattern, txPattern, identifier)
    if feature == 'gene':
        return DictParser(mode='bed6', feature=feature, save=save, **geneDict)
    else:
        return DictParser(mode='bed6', feature=feature, save=save, **txDict)


def FormatParser(attribute, formatType):
    formatDict = defaultdict(dict)
    listA = []
    if formatType == 'gff3':
        listA = re.split('=|;', attribute)
    elif formatType == 'gtf':
        listA = re.split('\s"|";\s', attribute)
    return MultiSys.List2Dict(listA)


def MainParser(Object, formatType, genePattern, txPattern, identifier):
    genePatternList = genePattern.split(',')
    txPatternList = txPattern.split(',')
    parserGeneDict = defaultdict(dict)
    parserTxDict = defaultdict(dict)
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
    annoRegex = re.compile(r'(^#.*)|(^$)')
    for line in lineList:
        if annoRegex.search(line) is not None:
            continue
        contList = line.rstrip('\n').split('\t')
        feature = contList[2]  # gene, transcript, exon, CDS...
        if feature == 'chromosome':
            continue
        chrom = contList[0]
        chrStart = int(contList[3]) - 1
        chrEnd = int(contList[4])
        strand = contList[6]
        attributeDict = FormatParser(contList[-1], formatType)
        parentID = idPattern = ''
        #
        if formatType == 'gff3':
            idPattern, parentPattern = identifier.split(',')
            if feature == 'CDS' or feature == 'exon':
                parentID = attributeDict[parentPattern]
            elif feature in txPatternList:
                parentID = attributeDict[parentPattern]
                idAttribute = attributeDict[idPattern]
            elif feature in genePatternList:
                idAttribute = attributeDict[idPattern]
            else:
                continue
        elif formatType == 'gtf':
            geneID, transcriptID = identifier.split(',')
            if feature == 'CDS' or feature == 'exon':
                parentID = attributeDict[transcriptID]
            elif feature in txPatternList:
                idAttribute = attributeDict[transcriptID]
                parentID = attributeDict[geneID]
            elif feature in genePatternList:
                idAttribute = attributeDict[geneID]
            else:
                continue
        #
        if feature in genePatternList:
            parserGeneDict[idAttribute] = [chrom, chrStart, chrEnd, strand]
        elif feature in txPatternList:
            parserTxDict[idAttribute]['tx'] = [
                chrom, chrStart, chrEnd, strand]
            parserTxDict[idAttribute]['pa'] = parentID
        else:
            # if feature not in parserTxDict.get(parentID, {}):
            if feature not in parserTxDict[parentID]:
                parserTxDict[parentID].setdefault(feature, {}).setdefault(
                    's', []).append(chrStart)
                parserTxDict[parentID][feature].setdefault(
                    'e', []).append(chrEnd)
            else:
                parserTxDict[parentID][feature]['s'].append(chrStart)
                parserTxDict[parentID][feature]['e'].append(chrEnd)
    return (parserGeneDict, parserTxDict)


def DictParser(mode, feature, save, **paserDict):
    bedDict = defaultdict(dict)
    if mode == 'bed12':
        for txID in paserDict:
            if 'tx' not in paserDict[txID]:
                continue
            locus = paserDict[txID]['tx']
            startList = sorted(paserDict[txID]['exon']['s'])
            endList = sorted(paserDict[txID]['exon']['e'])
            thickStart = thickEnd = locus[1]
            blockCount = len(startList)
            bedBlockStartCol = ','.join(
                MultiSys.List2Str([x - thickStart for x in startList]))  # MultiSys.List2Str(list(map(lambda x: x[0]-x[1], zip(endList, startList))))
            bedBlockLengthCol = ','.join(
                MultiSys.List2Str([endList[i] - startList[i] for i in range(blockCount)]))
            if 'CDS' in paserDict[txID]:
                thickStart = min(paserDict[txID]['CDS']['s'])
                thickEnd = max(paserDict[txID]['CDS']['e'])
            if save == 'list':
                bedDict[txID] = list([locus[0], locus[1], locus[2], txID, 0, locus[3], thickStart, thickEnd,
                                      255, blockCount, bedBlockLengthCol, bedBlockStartCol])
            else:
                bedDict[txID] = '\t'.join(MultiSys.List2Str([locus[0], locus[1], locus[2], txID, 0, locus[3], thickStart, thickEnd,
                                                             255, blockCount, bedBlockLengthCol, bedBlockStartCol]))
    elif mode == 'bed6':
        if feature == 'gene':
            for geneID in paserDict:
                locus = paserDict[geneID]
                if save == 'list':
                    bedDict[geneID] = list(
                        [locus[0], locus[1], locus[2], geneID, 0, locus[3]])
                else:
                    bedDict[geneID] = '\t'.join(
                        MultiSys.List2Str([locus[0], locus[1], locus[2], geneID, 0, locus[3]]))
        elif feature == 'transcript':
            for txID in paserDict:
                locus = paserDict[txID]['tx']
                if save == 'list':
                    bedDict[txID] = list(
                        [locus[0], locus[1], locus[2], txID, 0, locus[3]])
                else:
                    bedDict[txID] = '\t'.join(
                        MultiSys.List2Str([locus[0], locus[1], locus[2], txID, 0, locus[3]]))
        elif feature == 'exon':
            for txID in paserDict:
                if 'tx' not in paserDict[txID]:
                    continue
                locus = paserDict[txID]['tx']
                startList = sorted(paserDict[txID]['exon']['s'])
                endList = sorted(paserDict[txID]['exon']['e'])
                for i in range(len(startList)):
                    exonID = txID + ":exon:" + str(i + 1)
                    if save == 'list':
                        bedDict[exonID] = list(
                            [locus[0], startList[i], endList[i], exonID, 0, locus[3]])
                    else:
                        bedDict[exonID] = '\t'.join(MultiSys.List2Str(
                            [locus[0], startList[i], endList[i], exonID, 0, locus[3]]))
        elif feature == 'CDS':
            for txID in paserDict:
                if 'tx' not in paserDict[txID]:
                    continue
                locus = paserDict[txID]['tx']
                if 'CDS' not in paserDict[txID]:
                    continue
                startList = sorted(paserDict[txID]['CDS']['s'])
                endList = sorted(paserDict[txID]['CDS']['e'])
                for i in range(len(startList)):
                    cdsID = txID + ":CDS:" + str(i + 1)
                    if save == 'list':
                        bedDict[cdsID] = list(
                            [locus[0], startList[i], endList[i], cdsID, 0, locus[3]])
                    else:
                        bedDict[cdsID] = '\t'.join(
                            MultiSys.List2Str([locus[0], startList[i], endList[i], cdsID, 0, locus[3]]))
    return(bedDict)


def Test():
    file = sys.argv[1]
    bedDict = ToBed12(file, save='string')
    for x in range(bedList):
        print(bedDict[x])

if __name__ == '__main__':
    Test()
