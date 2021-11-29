#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
import operator
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
parser.add_argument('--grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('--grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('-i', '--input', action='store', type=str,
                    required=True,
                    help='bed output from ecConvergeAaToSegBed.py')
parser.add_argument('-o', '--output', action='store', type=str,
                    required=True,
                    help='output directory')
parser.add_argument('-p', '--prefix', action='store', type=str,
                    default='ecDNA',
                    help='prefix of output files')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

if bool(args.grepKept):
    kept = True
    regex = re.compile(r'{0}'.format(args.grepKept))
elif bool(args.grepExpel):
    kept = False
    regex = re.compile(r'{0}'.format(args.grepExpel))
else:
    regex = False

cycleClassDict = defaultdict(dict)
sampleCyDict = defaultdict(dict)
scInfoDict = defaultdict(dict)
with open(args.input, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        sampleName = row[6].replace('_encode', '')
        if bool(regex):
            if kept:
                if bool(regex.search(sampleName)) is False:
                    continue
            else:
                if bool(regex.search(sampleName)) is True:
                    continue
        cycleId = row[7]
        cycleLength = int(row[8])
        segNum = int(row[11])
        cycleClass = row[12]
        cycleClassDict['cycleClass'][cycleClass] = 1
        cycleClassDict['cycle'][cycleId] = cycleClass
        if sampleName not in sampleCyDict:
            sampleCyDict[sampleName] = defaultdict(set)
        sampleCyDict[sampleName][cycleClass].update([cycleId])
        if cycleId not in scInfoDict[sampleName]:
            scInfoDict[sampleName][cycleId] = defaultdict(dict)
            scInfoDict[sampleName][cycleId]['segNum'] = segNum
            scInfoDict[sampleName][cycleId]['length'] = cycleLength
            scInfoDict[sampleName][cycleId]['cover'] = defaultdict(list)
            scInfoDict[sampleName][cycleId]['geneClass'] = defaultdict(int)
            scInfoDict[sampleName][cycleId]['gene'] = defaultdict(int)

selectRepeatClassList = ["SINE", "LINE", "LTR", "Low_complexity", "Simple_repeat", "Satellite"]
geneClassList = list()
sampleGeneDict = defaultdict(dict)
elementList = list()
nonGeneList = ['repeat', 'microsatellite', 'CpG-island']
ecAnnoDict = defaultdict(list)
with open(args.anno, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        sampleName = row[6].replace('_encode', '')
        if bool(regex):
            if kept:
                if bool(regex.search(sampleName)) is False:
                    continue
            else:
                if bool(regex.search(sampleName)) is True:
                    continue
        cycleId = row[7]
        cycleClass = row[12]
        geneClass = row[18]
        geneType = row[19]
        geneName = row[15]
        chrom = row[0]
        if geneClass == 'intergenic':
            continue
        start = int(row[22])
        end = int(row[23])
        locus = [chrom, start, end]
        if geneClass in nonGeneList:
            if geneClass == 'repeat':
                element = geneName.split(':')[-1]
                if element == 'Alu':
                    geneClass = 'Alu'
                else:
                    if geneType in selectRepeatClassList:
                        geneClass = geneType
            element = geneClass
        else:
            element = 'gene'
            scInfoDict[sampleName][cycleId]['gene'][geneName] += 1
            if cycleClass not in sampleGeneDict[geneName]:
                sampleGeneDict[geneName][cycleClass] = defaultdict(int)
            sampleGeneDict[geneName][cycleClass][sampleName] += 1
        if element not in elementList:
            elementList.append(element)
        if geneClass not in geneClassList:
            geneClassList.append(geneClass)
        scInfoDict[sampleName][cycleId]['geneClass'][geneClass] += 1
        scInfoDict[sampleName][cycleId]['cover'][element].append(locus)

# store the statistics about total length of an element that cover a cycle
scStatsDict = defaultdict(dict)
scCovDict = defaultdict(dict)
for sampleName in sorted(scInfoDict.keys()):
    for cycleClass in sorted(sampleCyDict[sampleName]):
        cycleIdList = sorted(sampleCyDict[sampleName][cycleClass])
        ## sample-cycleClass
        scStatsDict[sampleName][cycleClass] = defaultdict(dict)
        scStatsDict[sampleName][cycleClass]['geneClass'] = defaultdict(int)
        scStatsDict[sampleName][cycleClass]['gene'] = defaultdict(int)
        scStatsDict[sampleName][cycleClass]['length'] = 0
        scStatsDict[sampleName][cycleClass]['segNum'] = 0
        scStatsDict[sampleName][cycleClass]['count'] = len(cycleIdList)
        for cycleId in cycleIdList:
            scCovDict[sampleName][cycleId] = defaultdict(dict)
            tempDict = scInfoDict[sampleName][cycleId]
            cycleLength = tempDict['length']
            ## sample-cycleClass
            for geneName in tempDict['gene'].keys():
                scStatsDict[sampleName][cycleClass]['gene'][geneName] += 1
            scStatsDict[sampleName][cycleClass]['length'] += cycleLength
            scStatsDict[sampleName][cycleClass]['segNum'] += tempDict['segNum']
            ## geneClass
            for geneClass in tempDict['geneClass'].keys():
                scStatsDict[sampleName][cycleClass]['geneClass'][geneClass] += tempDict['geneClass'][geneClass]
            ## element 
            for element in tempDict['cover'].keys():
                coverList = sorted(tempDict['cover'][element], key = operator.itemgetter(1, 2))
                newCoverList = list()
                plocus = coverList[0]
                coverListLen = len(coverList)
                if coverListLen > 1:
                    for i in range(1,coverListLen):
                        nlocus = coverList[i]
                        merge = bedutils.bedops(plocus, nlocus).merge()
                        if merge.mbool is False:
                            newCoverList.append(plocus)
                            plocus = nlocus
                            if (i +1) == coverListLen:
                                newCoverList.append(nlocus)
                        else:
                            plocus = [merge.m.chr, merge.m.start, merge.m.end]
                            if (i +1) == coverListLen:
                                newCoverList.append(plocus)
                else:
                    newCoverList.append(plocus)
                coverLength = sum(map(lambda x:x[2] - x[1], newCoverList))
                scCovDict[sampleName][cycleId][element] = coverLength

# output statistics
os.makedirs(args.output, exist_ok=True)

# sample-cycle
prefix = os.path.join(args.output, args.prefix + '.sample.')
cycleClassList = sorted(cycleClassDict['cycleClass'].keys())
typeList = ['count', 'length', 'segNum', 'gene']
for otype in typeList:
    outputFile = prefix + otype +'.stats'
    with open(outputFile, 'w') as output:
        row = ['sample']
        row.extend(cycleClassList)
        output.write('\t'.join(row) + '\n')
        for sampleName in sorted(scStatsDict.keys()):
            row = [sampleName]
            for cycleClass in cycleClassList:
                if cycleClass not in scStatsDict[sampleName]:
                    number = 0
                else:
                    tempDict = scStatsDict[sampleName][cycleClass]
                    number = tempDict[otype]
                    if otype != 'count':
                        if otype == 'gene':
                            number = len(tempDict['gene'].keys())
                        else:
                            count = tempDict['count']
                            number = round(number / count, 1)
                if otype == 'length':
                    number = number / 1000
                row.append(str(number))
            output.write('\t'.join(row) + '\n')

# store the statistics about 
geneClassList = sorted(geneClassList)
outputFile = prefix + 'geneClass.stats'
with open(outputFile, 'w') as output:
    row = ['sample', 'cycleClass']
    row.extend(geneClassList)
    output.write('\t'.join(row) + '\n')
    for sampleName in sorted(scStatsDict.keys()):
        reDict = defaultdict(dict)
        totalCount = 0
        totalNumList = [0 for i in range(len(geneClassList))]
        for cycleClass in cycleClassList:
            row = [sampleName, cycleClass]
            if cycleClass not in scStatsDict[sampleName]:
                for geneClass in geneClassList:
                    row.append('0')
            else:
                tempDict = scStatsDict[sampleName][cycleClass]
                count = tempDict['count']
                totalCount += count
                for i in range(len(geneClassList)):
                    geneClass = geneClassList[i]
                    if geneClass not in tempDict['geneClass']:
                        number = 0
                        pct = 0
                    else:
                        number = tempDict['geneClass'][geneClass]
                        pct = round(number / count, 1)
                    totalNumList[i] += number
                    row.append(str(pct))
                output.write('\t'.join(row) + '\n')


# sample-cycle-element coverage statistics
elementList = sorted(elementList)
cycleEleCovFile = os.path.join(args.output, args.prefix + '.sample.cycle.element.cov.stats')
sampleEleCovFile = os.path.join(args.output, args.prefix + '.sample.element.cov.stats')
with open(cycleEleCovFile, 'w') as cycleEleCov, open(sampleEleCovFile, 'w') as sampleEleCov:
    # for cycle-element
    cycleEleRow = ['sample', 'cycleId', 'cycleClass', 'cycleLength']
    cycleEleRow.extend(elementList)
    cycleEleCov.write('\t'.join(cycleEleRow) + '\n')
    # for sample-element
    sampleEleCovRow = ['sample', 'cycleClass']
    sampleEleCovRow.extend(elementList)
    sampleEleCov.write('\t'.join(sampleEleCovRow) + '\n')
    for sampleName in sorted(scCovDict.keys()):
        for cycleClass in sorted(sampleCyDict[sampleName].keys()):
            sampleEleCovRow = [sampleName, cycleClass]
            eleLengthList = [0] * len(elementList)
            totalCycleLength = 0
            cycleIdList= sorted(sampleCyDict[sampleName][cycleClass])
            for cycleId in cycleIdList:
                cycleLength = scInfoDict[sampleName][cycleId]['length']
                totalCycleLength += cycleLength
                cycleEleRow = [sampleName, cycleId, cycleClass, str(cycleLength)]
                for i in range(len(elementList)):
                    element = elementList[i]
                    if element not in scCovDict[sampleName][cycleId]:
                        coverLength = 0
                        coverage = 0
                    else:
                        coverLength = scCovDict[sampleName][cycleId][element]
                        coverage = round(coverLength / cycleLength * 100, 3)
                    cycleEleRow.append(str(coverage))
                    eleLengthList[i] += coverLength
                cycleEleCov.write('\t'.join(cycleEleRow) + '\n')
            for i in range(len(eleLengthList)):
                eleLengthList[i] = str(round( eleLengthList[i] / totalCycleLength * 100, 3))
            sampleEleCovRow += eleLengthList
            sampleEleCov.write('\t'.join(sampleEleCovRow) + '\n')

# gene-sample
geneSampleFile = os.path.join(args.output, args.prefix + '.gene.sample.stats')
with open(geneSampleFile, 'w') as output:
    row = ['gene', 'all']
    row.extend(cycleClassList)
    output.write('\t'.join(row) + '\n')
    for geneName in sorted(sampleGeneDict.keys()):
        row = [geneName]
        cycleClassRow = list()
        geneSampleList = list()
        for cycleClass in cycleClassList:
            if cycleClass not in sampleGeneDict[geneName]:
                sampleNum = '0'
            else:
                sampleList = sampleGeneDict[geneName][cycleClass]
                sampleNum = str(len(sampleList))
                geneSampleList.extend(sampleList)
            cycleClassRow.append(sampleNum)
        totalSampleNum = str(len(set(geneSampleList)))
        row.append(totalSampleNum)
        row.extend(cycleClassRow)
        output.write('\t'.join(row) + '\n')
