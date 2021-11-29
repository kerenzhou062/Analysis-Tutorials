#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
from pybedtools import BedTool
from copy import copy
import re
# in-house module
import bedutils

parser = argparse.ArgumentParser(
    description="This script is used for transforming coordinates from \
    .rebuild.txt to circle-based bed and fasta",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a', '--anno', action='store', type=str,
                    required=True,
                    help='*.rebuild.txt output from ecRebuildSegAnno.py')
parser.add_argument('-e', '--extend', action='store', type=int,
                    default=500,
                    help='get extended sequence around junction sites')
parser.add_argument('-i', '--input', action='store', type=str,
                    required=True,
                    help='bed output from ecConvergeAaToSegBed.py')
parser.add_argument('-f', '--fasta', action='store', type=str,
                    required=True,
                    help='genome fasta')
parser.add_argument('-o', '--output', action='store', type=str,
                    required=True,
                    help='output directory (prefix with sample name)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def chunkString(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

kwargs = {'name':True, 'tab':True, 's':True}

ecBedDict = defaultdict(dict)
with open(args.input, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        name = row[3]
        ecBedDict[name] = row

ecAnnoDict = defaultdict(list)
with open(args.anno, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        name = row[3]
        ecAnnoDict[name].append(row)

ecBed = BedTool(args.input)
ecSequence = ecBed.sequence(fi=args.fasta, **kwargs)
sampleFaDict = defaultdict(dict)
with open(ecSequence.seqfn, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        name = row[0].replace('(+)', '').replace('(-)', '')
        if name not in ecBedDict:
            continue
        sequence = row[1]
        nameRow = name.split('|')
        sampleName = nameRow[0]
        cycleName = '|'.join(nameRow[0:3])
        segOrder = int(nameRow[-1].split('-')[1])
        if sampleName not in sampleFaDict:
            sampleFaDict[sampleName] = defaultdict(dict)
        sampleFaDict[sampleName][cycleName][segOrder] = sequence

nonGeneClassList = ['repeat', 'microsatellite', 'CpG-island']
os.makedirs(args.output, exist_ok=True)
faDir = os.path.join(args.output, 'fa')
bedDir = os.path.join(args.output, 'bed')
os.makedirs(faDir, exist_ok=True)
os.makedirs(bedDir, exist_ok=True)

juncHeadRow = ['#rearrangement', 'start', 'end', 'junctionName', 'score', 'strand', 'rearrangementId', 'upstream', 'downstream']
geneBedHeadRow = ['#rearrangement', 'start', 'end', 'uniqName', 'score', 'strand', 'geneId', 'geneName', 'synonymous', 'description']
geneBedHeadRow += ['geneClass', 'geneType', 'transcriptId', 'transcriptType', 'rearrangementClass', 'segmentId']
additionBedHeadRow = ['#rearrangement', 'start', 'end', 'uniqName', 'score', 'strand', 'id', 'name', 'synonymous', 'description']
additionBedHeadRow += ['class', 'type', 'id', 'type', 'rearrangementClass', 'segmentId']

for sampleName in sorted(sampleFaDict.keys()):
    ecFasta = os.path.join(faDir, sampleName + '.fa')
    ecJuncFasta = os.path.join(bedDir, sampleName + '.junction.sequence.txt')
    ecJunctionBed = os.path.join(bedDir, sampleName + '.junction.bed')
    ecReBedGeneAnnoBed = os.path.join(bedDir, sampleName + '.anno.gene.bed')
    ecReBedStrandAnnoBed = os.path.join(bedDir, sampleName + '.anno.element.bed')
    ecReBedUnstrandAnnoBed = os.path.join(bedDir, sampleName + '.anno.unstrand.bed')
    with open(ecFasta, 'w') as fasta, open(ecJunctionBed, 'w') as juncBed, \
    open(ecJuncFasta, 'w') as juncSeqFasta, open(ecReBedGeneAnnoBed, 'w') as geneAnnoBed, \
    open(ecReBedStrandAnnoBed, 'w') as strandAnnoBed, open(ecReBedUnstrandAnnoBed, 'w') as unstrandAnnoBed:
        juncBed.write('\t'.join(juncHeadRow) + '\n')
        geneAnnoBed.write('\t'.join(geneBedHeadRow) + '\n')
        strandAnnoBed.write('\t'.join(additionBedHeadRow) + '\n')
        unstrandAnnoBed.write('\t'.join(additionBedHeadRow) + '\n')
        cycleCount = 1
        juncSeqList = list()
        for cycleName in sorted(sampleFaDict[sampleName].keys()):
            tempDict = sampleFaDict[sampleName][cycleName]
            orderList = sorted(tempDict.keys())
            sequenceList = [ tempDict[x] for x in orderList ]
            seqLength = 0
            newChrName = 'chr' + str(cycleCount)
            orderNum = len(orderList)
            for order in orderList:
                cycleSegId = '|'.join([cycleName, 'seg-' + str(order)])
                bedRow = copy(ecBedDict[cycleSegId])
                bedStart = int(bedRow[1])
                bedEnd = int(bedRow[2])
                segLength = int(bedRow[10])
                cycleClass = bedRow[-1]
                ## record junction sites and up-, down-stream sequences
                juncRowList = list()
                ## get junction sites
                juncStart = seqLength + segLength
                juncEnd = juncStart + 1
                ## get junction sites sequences
                if orderNum == 1:
                    upSeq = 'NA'
                    downSeq = 'NA'
                else:
                    if order == orderNum:
                        upSeq = 'NA'
                        downSeq = 'NA'
                    else:
                        seqStart = segLength - args.extend
                        upSeq = sequenceList[order - 1][seqStart:segLength]
                        downSeq = sequenceList[order][0:args.extend]
                juncName = '|'.join([cycleClass, newChrName, str(order), str(order + 1)])
                if orderNum != 1:
                    ## append row to junc
                    juncRow = [newChrName, str(juncStart), str(juncEnd), juncName, '255', '.', cycleName, upSeq, downSeq]
                    if upSeq != 'NA':
                        juncSeqList.append([juncName, upSeq, downSeq])
                    juncRowList.append(juncRow)
                if cycleClass == 'circular':
                    if order == 1:
                        ## junction sites sequence
                        juncSegLength = len(sequenceList[-1])
                        seqStart = juncSegLength - args.extend
                        upSeq = sequenceList[-1][seqStart:juncSegLength]
                        downSeq = sequenceList[0][0:args.extend]
                        ## append row to junc
                        juncName = '|'.join([cycleClass, newChrName, '0', '-1'])
                        juncSeqList.append([juncName, upSeq, downSeq])
                        juncRow = [newChrName, str(0), str(1), juncName, '255', '.', cycleName, upSeq, downSeq]
                        juncRowList.append(juncRow)
                    elif order == orderNum:
                        juncRowList = list()
                        juncStart = seqLength + segLength - 1
                        juncEnd = juncStart + 1
                        juncName = '|'.join([cycleClass, newChrName, '-1', '0'])
                        juncRow = [newChrName, str(juncStart), str(juncEnd), juncName, '255', '.', cycleName, 'NA', 'NA']
                        juncRowList.append(juncRow)
                for juncRow in juncRowList:
                    juncBed.write('\t'.join(juncRow) + '\n')
                ## re-annotate gene annotation
                geneCount = 1
                for annoRow in ecAnnoDict[cycleSegId]:
                    geneType = annoRow[19]
                    if geneType == 'intergenic':
                        continue
                    annoStart = int(annoRow[22])
                    annoEnd = int(annoRow[23])
                    annoLength = annoEnd - annoStart
                    annoInfoList = annoRow[14:22]
                    annoInfoList.append(cycleClass)
                    annoInfoList.append(cycleSegId)
                    orient = annoRow[27]
                    geneId = annoRow[14].split('.')[0]
                    geneName = annoRow[15]
                    bedName = '|'.join([geneName, newChrName, geneId, str(order), str(geneCount)])
                    geneClass = annoRow[18]
                    if bedStart <= annoStart and annoEnd <= bedEnd:
                        annoStart = seqLength + annoStart - bedStart
                        annoEnd = annoStart + annoLength
                        newAnnoRow = [newChrName, str(annoStart), str(annoEnd), bedName, '255', orient]
                        newAnnoRow.extend(annoInfoList)
                        if orient == '.':
                            unstrandAnnoBed.write('\t'.join(newAnnoRow) + '\n')
                        else:
                            if geneClass in nonGeneClassList:
                                strandAnnoBed.write('\t'.join(newAnnoRow) + '\n')
                            else:
                                geneAnnoBed.write('\t'.join(newAnnoRow) + '\n')
                    geneCount += 1
                seqLength += segLength
            ## ecDNA sequence
            sequence = ''.join(sequenceList)
            sequenceList = list(chunkString(sequence, 50))
            sequenceList[-1] += '\n'
            fasta.write('\n'.join(['>' + newChrName, '\n'.join(sequenceList)]))
            cycleCount += 1
        ## ecDNA junction sites
        juncSeqList = sorted(juncSeqList, key=lambda x:x[0])
        for juncSeqRow in juncSeqList:
            upSeq = '\n'.join(list(chunkString(juncSeqRow[1].upper(), 50)))
            downSeq = '\n'.join(list(chunkString(juncSeqRow[2].upper(), 50)))
            juncSeqFasta.write('>{0}|up\n{1}\n'.format(juncSeqRow[0], upSeq))
            juncSeqFasta.write('>{0}|down\n{1}\n'.format(juncSeqRow[0], downSeq))
