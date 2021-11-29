#!/usr/bin/env python3
import os
import sys
import argparse
from glob import glob
from copy import copy
import re
from collections import defaultdict
import pybedtools

parser = argparse.ArgumentParser(
    description="This script is used for converging ecDNA results from runPipeAA.sh to bed6+",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-c', '--ecdna', action='store', type=str, required=True,
                    help='ecDNA bed')
parser.add_argument('-e', '--ecpeak', action='store', type=str, required=True,
                    help='files that ecDNA bed intersected with ChIP-seq peaks')
parser.add_argument('-r', '--repeak', action='store', type=str, required=True,
                    help='files that repeat bed intersected with ChIP-seq peaks')
parser.add_argument('-s', '--sample', action='store', type=str, required=True,
                    help='sample name of input --ecpeak')
parser.add_argument('-p', '--prefix', action='store', type=str, required=True,
                    help='prefix for output files')
parser.add_argument('--cut1', action='store', type=float,
                    default=15,
                    help='cutoff for filtering (eg. ecDNA-like >= cut1)')
parser.add_argument('--cut2', action='store', type=float,
                    default=2,
                    help='cutoff for filtering (ecDNA-like >= cut2 * linear & >= cut2 * No_fSCNA)')
parser.add_argument('--detail', action='store_true',
                    default=False,
                    help='run with detail mode of repeat element')
parser.add_argument('-o', '--output', action='store', type=str, required=True,
                    help='output directory')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

sampleDict = defaultdict(int)
with open(args.ecdna, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        sample = row[6]
        sampleDict[sample] += 1

if args.sample not in sampleDict:
    print("--sample is not valid!")
    print("Available samples are listed as follow:")
    for sample in sorted(sampleDict.keys()):
        print(sample)
    sys.exit(0)

if args.detail is False:
    selectRepeatClassList = ["Alu", "SINE", "LINE", "LTR", "Low_complexity", "Simple_repeat", "Satellite"]
else:
    selectRepeatClassList = set()
factorPeakReaptDict = defaultdict(dict)
with open(args.repeak, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        repeatRow = row[0:6]
        repeatId = repeatRow[3]
        if bool(re.search(r'\?', repeatId)):
            continue
        repeatIdRow = repeatId.split(':')
        if args.detail is False:
            if repeatIdRow[2] == "Alu":
                repeatClass = "Alu"
            else:
                repeatClass = repeatIdRow[1]
            ## only keep peaks overlap with Alu or SINE or LINE element
            if repeatClass not in selectRepeatClassList:
                continue
        else:
            if repeatIdRow[1] == repeatIdRow[2]:
                repeatClass = repeatRow[2]
            else:
                repeatIdRow[2] = repeatIdRow[2].replace('-', '_')
                repeatClass = '_'.join(repeatIdRow[1:3])
                selectRepeatClassList.update([repeatClass])
        peakRow = row[6:]
        peakId = peakRow[3]
        factor = peakRow[-1]
        # record factor -> rpeatid
        if peakId not in factorPeakReaptDict[factor]:
            factorPeakReaptDict[factor][peakId] = defaultdict(int)
        factorPeakReaptDict[factor][peakId][repeatClass] += 1

selectRepeatClassList = sorted(selectRepeatClassList)

ecdnaDict = defaultdict(dict)
ecdnaDict['seg'] = defaultdict(dict)
ecdnaDict['cycle'] = defaultdict(dict)
ecdnaDict['sample'] = defaultdict(set)
with open(args.ecdna, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        segId = row[3]
        sample = row[6]
        if sample != args.sample:
            continue
        cycleId = row[7]
        ## cycle length in bp
        cycleLength = int(row[8])
        cycleClass = row[12]
        cycleClassCustom = row[12]
        ecdnaDict['seg'][segId] = cycleId
        if cycleId not in ecdnaDict['cycle']:
            ecdnaDict['cycle'][cycleId]['seg'] = list()
        ecdnaDict['cycle'][cycleId]['seg'].append(segId)
        ecdnaDict['cycle'][cycleId]['cycleClass'] = cycleClass
        ecdnaDict['cycle'][cycleId]['cycleLength'] = cycleLength
        ecdnaDict['cycle'][cycleId]['cycleClassCustom'] = cycleClassCustom
        ecdnaDict['sample'][cycleClass].update([cycleId])

segIdLinkDict = defaultdict(str)
factorPeakDict = defaultdict(dict)
segIdFactorDict = defaultdict(dict)
with open(args.ecpeak, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        # fragment information
        # chr19, 21775969, 21786507, segment-44, 2, -, 21780969, 21781507, K562_encode|amplicon144|cycle-10|seg-6,K562_encode|amplicon144|cycle-6|seg-3,
        # Rearranged,Rearranged, translocation,translocation
        ecdnaSegRow = row[0:11]
        segReformId = ecdnaSegRow[3]
        ecStart = int(ecdnaSegRow[6])
        ecEnd = int(ecdnaSegRow[7])
        uniqSegIdList = ecdnaSegRow[8].split(',')
        for uniqSegId in uniqSegIdList:
            segIdLinkDict[uniqSegId] = segReformId
        # peak information
        peakRow = row[11:]
        peakId = peakRow[3]
        peakStart = int(peakRow[1])
        peakEnd = int(peakRow[2])
        factor = peakRow[-1]
        ## only consider peaks located in the rearranged fragment
        if ecStart <= peakStart and peakEnd <= ecEnd:
            if peakId not in factorPeakDict[factor]:
                factorPeakDict[factor][peakId] = '\t'.join(peakRow)
            if segReformId not in segIdFactorDict:
                segIdFactorDict[segReformId] = defaultdict(list)
            segIdFactorDict[segReformId][factor].append(peakId)

# build dictionary: cycle -> factor -> [peakIdList]
cycleFactorDict = defaultdict(dict)
for cycleClass in sorted(ecdnaDict['sample'].keys()):
    for cycleId in ecdnaDict['sample'][cycleClass]:
        if cycleId not in cycleFactorDict:
            cycleFactorDict[cycleId] = defaultdict(list)
        segReformIdList = list()
        for segId in sorted(ecdnaDict['cycle'][cycleId]['seg']):
            if segId in segIdLinkDict:
                segReformIdList.append(segIdLinkDict[segId])
        if len(segReformIdList) == 0:
            continue
        for segReformId in segReformIdList:
            for factor in sorted(segIdFactorDict[segReformId].keys()):
                cycleFactorDict[cycleId][factor].extend(segIdFactorDict[segReformId][factor])

# store the statistics about the peak number of a factor that binds to a cycle and repeat elements
cycleFactorStatsFile = os.path.join(args.output, args.prefix + '.cycle_factor_stats.txt')
factorCycleStatsDict = defaultdict(dict)
cycleClassStatsDict = defaultdict(dict)
with open(cycleFactorStatsFile, 'w') as out:
    row = ["cycleId", "cycleClass", "cycleLength", "factor", "peakIdNum", "peakIdDensity", "peakInRepeatNum", "repeatNum"] + selectRepeatClassList
    out.write('\t'.join(map(str, row)) + '\n')
    for cycleClass in sorted(ecdnaDict['sample'].keys()):
        for cycleId in sorted(ecdnaDict['sample'][cycleClass]):
            cycleLength = ecdnaDict['cycle'][cycleId]['cycleLength']
            if cycleId not in cycleFactorDict:
                continue
            for factor in sorted(cycleFactorDict[cycleId].keys(), key=lambda x:len(cycleFactorDict[cycleId][x]), reverse=True):
                # peaks that bind to fragment
                peakIdList = cycleFactorDict[cycleId][factor]
                peakIdNum = len(peakIdList)
                # peak number per kb
                peakIdDensity = peakIdNum / cycleLength * 1000
                # number of peaks that overlap with repeat regions
                repeatIdList = []
                peakInRepeatNum = 0
                for peakId in peakIdList:
                    if peakId in factorPeakReaptDict[factor]:
                        peakInRepeatNum += 1
                        repeatIdList.extend(factorPeakReaptDict[factor][peakId])
                repeatNum = len(repeatIdList)
                # peaks that overlap with a specific type of repeat regions
                peakInRepeatNumList = [0] * len(selectRepeatClassList)
                peakInRepeatPercentageList = []
                for i in range(len(selectRepeatClassList)):
                    repeatClass = selectRepeatClassList[i]
                    for peakId in peakIdList:
                        if peakId in factorPeakReaptDict[factor]:
                            if repeatClass in factorPeakReaptDict[factor][peakId]:
                                peakInRepeatNumList[i] += 1
                # percentage of peaks that overlap with a specific type of repeat regions
                peakInRepeatPercentageList = list(map(lambda x: x/peakIdNum * 100, peakInRepeatNumList))
                # output row
                row = [cycleId, cycleClass, cycleLength, factor, len(peakIdList), peakIdDensity, peakInRepeatNum, repeatNum] + peakInRepeatPercentageList
                out.write('\t'.join(map(str, row)) + '\n')
                # store the statistic information
                if cycleClass not in factorCycleStatsDict[factor]:
                    factorCycleStatsDict[factor][cycleClass] = defaultdict(dict)
                    factorCycleStatsDict[factor][cycleClass]['cycleId'] = []
                    factorCycleStatsDict[factor][cycleClass]['count'] = 0
                    factorCycleStatsDict[factor][cycleClass]['peakNum'] = 0
                    factorCycleStatsDict[factor][cycleClass]['peakInRepeatNum'] = [0] * len(selectRepeatClassList)
                factorCycleStatsDict[factor][cycleClass]['cycleId'].append(cycleId)
                factorCycleStatsDict[factor][cycleClass]['count'] += 1
                factorCycleStatsDict[factor][cycleClass]['peakNum'] += peakIdNum
                for i in range(len(selectRepeatClassList)):
                    factorCycleStatsDict[factor][cycleClass]['peakInRepeatNum'][i] += peakInRepeatNumList[i]


# store the statistics about how many times that a factor occur in a cycleClass, the peak number of a factor that binds to cycleClass and repeat elements
factorCycleRepeatClassStatsFile = os.path.join(args.output, args.prefix + '.factor_cycleClass_repeat_stats.txt')
factorEleCycleFilterFile = os.path.join(args.output, args.prefix + '.factor_repeat_filter.txt')
cycleClassList = sorted(ecdnaDict['sample'].keys())
with open(factorCycleRepeatClassStatsFile, 'w') as statsOut, open(factorEleCycleFilterFile, 'w') as filterOut:
    statsRow = ['factor', 'cycleClass', 'cycleCount', 'factorInCycleCount', 'factorIncyclePercentage', 'totalPeakNum', 'peakDensity'] + selectRepeatClassList
    statsOut.write('\t'.join(statsRow) + '\n')
    filterRow = ['repeat'] + statsRow
    filterOut.write('\t'.join(filterRow) + '\n')
    for factor in sorted(factorCycleStatsDict.keys()):
        eleCyclePercentageDict = defaultdict(dict)
        filterRowList = []
        for cycleClass in cycleClassList:
            totalCycleNum = len(ecdnaDict['sample'][cycleClass])
            statsRow = [factor, cycleClass]
            if cycleClass in factorCycleStatsDict[factor]:
                totalCycleLength = sum(map(lambda x: ecdnaDict['cycle'][x]['cycleLength'], set(factorCycleStatsDict[factor][cycleClass]['cycleId'])))
                factorInCycleCount = factorCycleStatsDict[factor][cycleClass]['count']
                factorIncyclePercentage = factorInCycleCount / totalCycleNum
                totalPeakNum = factorCycleStatsDict[factor][cycleClass]['peakNum']
                # peak density per kb
                peakDensity = totalPeakNum / totalCycleLength * 1000 * factorIncyclePercentage
                factorIncyclePercentage = factorIncyclePercentage * 100
                peakInRepeatNumList = factorCycleStatsDict[factor][cycleClass]['peakInRepeatNum']
                peakInRepeatPercentageList = list(map(lambda x: str(x/totalPeakNum * 100), peakInRepeatNumList))
            else:
                factorInCycleCount = 0
                factorIncyclePercentage = 0
                peakDensity = 0
                peakInRepeatPercentageList = ['0'] * len(selectRepeatClassList)
            for i in range(len(selectRepeatClassList)):
                repeatClass = selectRepeatClassList[i]
                peakInRepeatPercentage = peakInRepeatPercentageList[i]
                if cycleClass not in eleCyclePercentageDict[repeatClass]:
                    eleCyclePercentageDict[repeatClass][cycleClass] = defaultdict(dict)
                eleCyclePercentageDict[repeatClass][cycleClass]['totalPeakNum'] = totalPeakNum
                eleCyclePercentageDict[repeatClass][cycleClass]['peakInCyclePc'] = factorIncyclePercentage
                eleCyclePercentageDict[repeatClass][cycleClass]['peakRepeatPc'] = float(peakInRepeatPercentage)
            # output row
            statsRow.extend(list(map(str, [totalCycleNum, factorInCycleCount, factorIncyclePercentage, totalPeakNum, peakDensity])))
            statsRow += peakInRepeatPercentageList
            statsOut.write('\t'.join(statsRow) + '\n')
            filterRowList.append(statsRow)
        # output row
        for repeatClass in sorted(eleCyclePercentageDict.keys()):
            ecdnaPeakInCyclePc = eleCyclePercentageDict[repeatClass]['ecDNA-like']['peakInCyclePc']
            ecdnaPeakTotalNum = eleCyclePercentageDict[repeatClass]['ecDNA-like']['totalPeakNum']
            ecdnaCycleNum = len(ecdnaDict['sample'][cycleClass])
            ecdnaPeakRepeatPc = eleCyclePercentageDict[repeatClass]['ecDNA-like']['peakRepeatPc']
            linearPeakRepeatPc = eleCyclePercentageDict[repeatClass]['Linear']['peakRepeatPc']
            nofscdnaPeakRepeatPc = eleCyclePercentageDict[repeatClass]['No_fSCNA']['peakRepeatPc']
            # filter
            if ecdnaPeakInCyclePc < 25 or ecdnaPeakTotalNum < ecdnaCycleNum:
                continue
            if ecdnaPeakRepeatPc < args.cut1:
                continue
            skipFlag = False
            for repeatClassCandidate in sorted(eleCyclePercentageDict.keys()):
                if repeatClass == repeatClassCandidate:
                    continue
                if eleCyclePercentageDict[repeatClass]['ecDNA-like']['peakRepeatPc'] < eleCyclePercentageDict[repeatClassCandidate]['ecDNA-like']['peakRepeatPc']:
                    skipFlag = True
                    break
            if skipFlag is True:
                continue
            if ecdnaPeakRepeatPc >= args.cut2 * linearPeakRepeatPc and ecdnaPeakRepeatPc >= args.cut2 * nofscdnaPeakRepeatPc:
                for filterRow in filterRowList:
                    filterRow = [repeatClass] + filterRow
                    filterOut.write('\t'.join(filterRow) + '\n')
