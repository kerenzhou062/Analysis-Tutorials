#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
from pybedtools import BedTool
# in-house module
import bedutils

parser = argparse.ArgumentParser(
    description="This script is used for annotating bed output from convergeAaToSegbedinfo.py",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a', '--anno', action='store', type=str,
                    required=True,
                    help='The annotation bed12 file')
parser.add_argument('-b', '--bin', action='store', type=str,
                    choices=['average', 'constant'],
                    default='constant',
                    help='The type of bin size')
parser.add_argument('-e', '--extend', action='store', type=int,
                    default=0,
                    help='Extend --extend bp around peak center ')
parser.add_argument('-i', '--input', nargs='+', type=str,
                    required=True,
                    help='The input files (bed3, bed6, bed12)')
parser.add_argument('-f', '--feature', nargs='+', type=str,
                    choices=["coding", "exon", "intron", "full"],
                    default='coding',
                    help='The bin features [coding:utr5,cds,utr3]')
parser.add_argument('-g', '--gene', nargs='+', type=str,
                    choices=["protein_coding", "non_coding", "all"],
                    default='protein_coding',
                    help='The bin features')
parser.add_argument('-m', '--method', action='store', type=str,
                    choices=['center', 'interval', 'exon'],
                    default='center',
                    help='The method for calculating bin coverage')
parser.add_argument('-n', '--name', nargs='+', type=str,
                    help='The sample name for --input files')
parser.add_argument('-o', '--output', action='store', type=str,
                    default="bin.txt",
                    help='The output file')
parser.add_argument('-p', '--peaknum', action='store', type=int,
                    help='Set the total number of peaks for --input')
parser.add_argument('-s', '--smooth', action='store', type=str,
                    choices=['none', 'move', 'average'],
                    default='move',
                    help='The smooth method')
parser.add_argument('-t', '--type', action='store', type=str,
                    choices=['density', 'percentage', 'number'],
                    default='density',
                    help='The type of output values')
parser.add_argument('-w', '--width', action='store', type=int,
                    default=5,
                    help='The span for smooth')
parser.add_argument('-z', '--size', action='store', type=int,
                    default=100,
                    help='The bin size of each feature')
parser.add_argument('--matchid', action='store_true',
                    default=False,
                    help='Match input bed name in annotation bed12')
parser.add_argument('--strand', action='store_true',
                    default=False,
                    help='Only report overlap on the same strand')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def ExonLenSum(blist):
    return sum(map(lambda x: x[1] - x[0], blist))

def RebuildBed(bedFile, method, extend):
    bedLineRow = list()
    regex = re.compile(r'^#')
    colNum = 0
    lineNum = 1
    with open(bedFile, 'r') as f:
        for line in f:
            if bool(regex.match(line)):
                continue
            else:
                row = line.strip().split('\t')
                bedinfo = bedutils.buildbed(row)
                uniqPeakName = '|'.join([bedinfo.name, str(lineNum)])
                if colNum == 0:
                    colNum = len(row)
                if method == 'center':
                    center = int((bedinfo.end + bedinfo.start) / 2)
                    row[1] = str(center - extend)
                    row[2] = str(center + extend + 1)
                    lineRow = row[0:3]
                    name = '##'.join([uniqPeakName, '1'])
                    lineRow.extend([name, str(bedinfo.score), bedinfo.strand])
                    bedLineRow.append('\t'.join(lineRow))
                elif method == 'exon':
                    bed12 = bedinfo.decode()
                    score = bed12.exon
                    for i in range(len(bed12.exon)):
                        exonId = str(i + 1)
                        name = '##'.join([uniqPeakName, exonId])
                        lineRow = [bed12.chr, bed12.exon[i][0], bed12.exon[i][1]]
                        lineRow.extend([name, bedinfo.score, bedinfo.strand])
                        bedLineRow.append('\t'.join(map(str,lineRow)))
                else:
                    if extend > 0:
                      center = int((bedinfo.end + bedinfo.start) / 2)
                      row[1] = center - extend
                      if row[1] < bedinfo.start:
                          row[1] = bedinfo.start
                      elif row[1] < 0:
                          row[1] = 0
                      else:
                          row[1] = center - extend
                      row[2] = center + extend + 1
                      if row[2] > bedinfo.end:
                          row[2] = bedinfo.end
                      row[1] = str(row[1])
                      row[2] = str(row[2])
                    lineRow = row[0:3]
                    name = '##'.join([uniqPeakName, '1'])
                    lineRow.extend([name, str(bedinfo.score), bedinfo.strand])
                    bedLineRow.append('\t'.join(lineRow))
                lineNum += 1
    bedLineRow.append('')
    bedlines = '\n'.join(bedLineRow)
    return [bedlines, lineNum]

def AnnoBed12ToBed6(bed12File, geneType, feature, binType, binsize):
    ## sub main
    bed12Dict = defaultdict(dict)
    statsDict = defaultdict(int)
    bed12InfoList = list()
    regex = re.compile(r'^#')
    count = 0
    with open(bed12File, 'r') as f:
        for line in f:
            if bool(regex.match(line)):
                continue
            row = line.strip().split('\t')
            bedinfo = bedutils.buildbed(row)
            bed12 = bedinfo.decode()
            bed12Dict[bedinfo.name]['info'] = bedinfo
            bed12Dict[bedinfo.name]['full'] = [[bedinfo.start, bedinfo.end]]
            bed12Dict[bedinfo.name]['exon'] = bed12.exon
            bed12Dict[bedinfo.name]['intron'] = bed12.intron
            ## length and size
            bed12Dict[bedinfo.name]['full_length'] = bed12.end - bed12.start
            bed12Dict[bedinfo.name]['intron_size'] = ExonLenSum(bed12.intron)
            bed12Dict[bedinfo.name]['exon_size'] = ExonLenSum(bed12.exon)

            if geneType == 'protein_coding':
                if len(bed12.cds) == 0:
                    continue
                bed12Dict[bedinfo.name]['utr5'] = bed12.utr5
                bed12Dict[bedinfo.name]['cds'] = bed12.cds
                bed12Dict[bedinfo.name]['utr3'] = bed12.utr3
            else:
                if len(bed12.cds) > 0:
                    continue
            ## calculate total length of features
            if binType == 'average' and feature == 'coding':
                for key in ['utr5', 'cds', 'utr3']:
                    statsDict[key] += ExonLenSum(bed12Dict[bedinfo.name][key])
            count += 1
    ## statistics, average lenth of features per transcript
    ## determin bin size of features
    binsizeDict = defaultdict(int)
    if feature == 'coding':
        featureNum = 3
    else:
        featureNum = 1
    binsizeTotal = featureNum * binsize
    if binType == 'average' and feature == 'coding':
        for key in sorted(statsDict.keys()):
            statsDict[key] = int(statsDict[key] / count)
        totalAveLen = sum([statsDict['utr5'], statsDict['cds'], statsDict['utr3']])
        utr5BinLen = round(statsDict['utr5'] / totalAveLen * binsizeTotal)
        cdsBinLen = round(statsDict['cds'] / totalAveLen * binsizeTotal)
        utr3BinLen = round(statsDict['utr3'] / totalAveLen * binsizeTotal)
        sumBinLen = sum([utr5BinLen, cdsBinLen, utr3BinLen])
        if sumBinLen < binsizeTotal:
            cdsBinLen += binsizeTotal - sumBinLen
        binsizeDict['utr5'] = utr5BinLen
        binsizeDict['cds'] = cdsBinLen
        binsizeDict['utr3'] = utr3BinLen
    else:
        if feature == 'coding':
            binsizeDict['utr5'] = binsize
            binsizeDict['cds'] = binsize
            binsizeDict['utr3'] = binsize
        else:
            binsizeDict[feature] = binsize
    ## create new annotation bed6
    ## fLenPreSum: comulative length of feature before this feature block
    bed6Dict = defaultdict(list)
    bedLineRow = list()
    for name in sorted(bed12Dict.keys()):
        tempDict = bed12Dict[name]
        bedinfo = tempDict['info']
        for feature in sorted(binsizeDict.keys()):
            if geneType == 'protein_coding':
                if 'cds' not in tempDict:
                    continue
            if len(tempDict[feature]) == 0:
                continue
            fLenTotal = ExonLenSum(tempDict[feature])
            mapLenPerbp = binsizeDict[feature] / fLenTotal
            fLenPreSum = 0
            for coord in tempDict[feature]:
                coordName = '##'.join([name, feature, str(fLenPreSum)])
                row = [bedinfo.chr, coord[0], coord[1], coordName, mapLenPerbp, bedinfo.strand]
                bed6Dict[coordName] = row
                bedLineRow.append('\t'.join(map(str, row)))
                fLenPreSum += coord[1] - coord[0]
    bedLineRow.append('')
    bedlines = '\n'.join(bedLineRow)
    return [bedlines, bed12Dict, bed6Dict, binsizeDict]

def Smooth(valueList, method, span):
    def span_ave(vlist, rstart, rend, rlen):
        return sum(map(lambda x:vlist[x], range(rstart, rend))) / rlen
    ## sub main
    vlength = len(valueList)
    smoothValList = list()
    if method == 'none':
        smoothValList = valueList
    elif method == 'move':
        for i in range(vlength):
            if i < span:
                rstart, rend, rlen= [0, i + 1, i + 1]
            else:
                rstart, rend, rlen= [i-span, i + 1, span]
            valueAve = span_ave(valueList, rstart, rend, rlen)
            smoothValList.append(valueAve)
    elif method == 'average':
        hspan = int(span / 2)
        for i in range(vlength):
            if i < hspan:
                tspan = i * 2 + 1
                rstart, rend, rlen= [0, i * 2 + 1, tspan]
            elif i >= (vlength - hspan):
                iToEnd = vlength - i
                tspan = (iToEnd * 2 -1) if (iToEnd % 2 == 1) else (iToEnd * 2)
                rstart, rend, rlen= [vlength - tspan, vlength, tspan]
            else:
                rstart, rend, rlen= [i - hspan, i + hspan + 1, span]
            valueAve = span_ave(valueList, rstart, rend, rlen)
            smoothValList.append(valueAve)
    return smoothValList

def RunMetagene(inputFile, annoBedDict, args, kwargs):
    annoBed = annoBedDict['bedtool']
    bed12Dict = annoBedDict['bed12']
    bed6Dict = annoBedDict['bed6']
    binsizeDict = annoBedDict['binsize']
    ## rebuild input bed
    inputBedLines, totalPeakNum = RebuildBed(inputFile, args.method, args.extend)
    if bool(args.peaknum):
        totalPeakNum = args.peaknum
    ## create bedtool objects
    inputBed = BedTool(inputBedLines, from_string=True).sort()
    intersect = inputBed.intersect(annoBed, **kwargs)
    interDict = defaultdict(dict)
    interStatsDict = defaultdict(dict)
    for item in intersect:
        overlapRow = item.fields[0:6]
        annoRow = item.fields[6:]
        peakName, childId = overlapRow[3].split('##')
        uniqFeatureName = annoRow[3]
        annoName = annoRow[3].split('##')[0]
        ## if set --matchid, force annoName containing peakName
        if args.matchid:
            originalPeakName = peakName.split("|")[0]
            if originalPeakName not in annoName:
                continue
        ## decode bed
        interStart = int(overlapRow[1])
        interEnd = int(overlapRow[2])
        interLen = interEnd - interStart
        if peakName not in interDict:
            interDict[peakName] = defaultdict(dict)
            interStatsDict[peakName] = defaultdict(int)
            if annoName not in interDict[peakName]:
                interDict[peakName][annoName] = defaultdict(dict)
        else:
            if annoName not in interDict[peakName]:
                interDict[peakName][annoName] = defaultdict(dict)
        interDict[peakName][annoName][childId] = [uniqFeatureName, interStart, interEnd]
        interStatsDict[peakName][annoName] += interLen
    
    ## determin unique peakName-annoName relationship
    ## intersection lengh -> longest transcript
    peakAnnoPairDict = defaultdict(str)
    for peakName in interStatsDict.keys():
        for annoName in interStatsDict[peakName].keys():
            if peakName not in peakAnnoPairDict:
                peakAnnoPairDict[peakName] = annoName
            else:
                preAnnoName = peakAnnoPairDict[peakName]
                preInterLen = interStatsDict[peakName][preAnnoName]
                interLen = interStatsDict[peakName][annoName]
                if interLen > preInterLen:
                    peakAnnoPairDict[peakName] = annoName
                elif interLen == preInterLen:
                    sizeKey = args.feature + '_size'
                    if args.feature == 'coding':
                        sizeKey = 'exon_size'
                    preAnnoLen = bed12Dict[preAnnoName][sizeKey]
                    annoLen = bed12Dict[annoName][sizeKey]
                    if annoLen > preAnnoLen:
                        peakAnnoPairDict[peakName] = annoName
    
    ## determin bin value
    binValDict = defaultdict(dict)
    for feature in binsizeDict.keys():
        binsize = binsizeDict[feature]
        binValDict[feature] = defaultdict(dict)
        for i in range(binsize):
            binValDict[feature][i]['sum'] = 0
            binValDict[feature][i]['peak'] = defaultdict(int)
    interTxDict = defaultdict(int)
    for peakName in sorted(peakAnnoPairDict.keys()):
        annoName = peakAnnoPairDict[peakName]
        for childId in interDict[peakName][annoName].keys():
            uniqFeatureName, interStart, interEnd = interDict[peakName][annoName][childId]
            interLen = interEnd - interStart
            annoName, feature, fLenPreSum = uniqFeatureName.split('##')
            fLenPreSum = int(fLenPreSum)
            annoStart = bed6Dict[uniqFeatureName][1]
            annoEnd = bed6Dict[uniqFeatureName][2]
            annoStrand = bed6Dict[uniqFeatureName][5]
            mapLenPerbp = bed6Dict[uniqFeatureName][4]
            interTxDict[annoName] += 1
            if annoStrand == '+':
                offset = interStart - annoStart + fLenPreSum
                for i in range(0,interLen):
                    binCoord = int((offset + i) * mapLenPerbp)
                    if binCoord >= binsizeDict[feature]:
                        binCoord = binsizeDict[feature] - 1
                    binValDict[feature][binCoord]['sum'] += 1
                    binValDict[feature][binCoord]['peak'][peakName] += 1
            else:
                offset = annoEnd - interEnd + fLenPreSum
                for i in range(-interLen, 0):
                    binCoord = int((offset - i) * mapLenPerbp)
                    if binCoord >= binsizeDict[feature]:
                        binCoord = binsizeDict[feature] - 1
                    binValDict[feature][binCoord]['sum'] += 1
                    binValDict[feature][binCoord]['peak'][peakName] += 1
    
    ## generate final bin-value dict
    interTxNum = len(interTxDict.keys())
    finalBinValDict = defaultdict(dict)
    featureList = [args.feature]
    count = 1
    if args.feature == 'coding':
        featureList = ['utr5', 'cds', 'utr3']
    for feature in featureList:
        binValList = list()
        for binCoord in sorted(binValDict[feature].keys()):
            binSum = binValDict[feature][binCoord]['sum']
            if args.type == 'density':
                binPeakNum = len(binValDict[feature][binCoord]['peak'].keys())
                binVal = binPeakNum / totalPeakNum * 100
            elif args.type == 'percentage':
                binVal = binSum / totalPeakNum * 100
            else:
                binVal = binSum
            binValList.append(binVal)
        smmothBinValList = Smooth(binValList, args.smooth, args.width)
        for binVal in smmothBinValList:
            if args.type != 'none':
                binVal = round(binVal, 5)
            finalBinValDict[count]['feature'] = feature
            finalBinValDict[count]['value'] = binVal
            count += 1
    return finalBinValDict

# main program
if __name__ == '__main__':
    kwargs = {'wb':True, "sorted":True, 's':args.strand}
    annoBedLines, bed12Dict, bed6Dict, binsizeDict = AnnoBed12ToBed6(args.anno, args.gene, args.feature, args.bin, args.size)
    annoBed = BedTool(annoBedLines, from_string=True).sort()
    annoBedDict = {'bedtool':annoBed, 'bed12':bed12Dict, 'bed6':bed6Dict, 'binsize':binsizeDict}
    sampleNameList = list()
    if bool(args.name):
        if len(args.input) != len(args.name):
            sys.error.write('--input and --name should match the number!')
            sys.exit()
        else:
            sampleNameList = args.name
    else:
        for inputFile in args.input:
            basename = os.path.splitext(os.path.split(inputFile)[-1])[0]
            sampleNameList.append(basename)
    binSampleValDict = defaultdict(dict)
    for i in range(len(args.input)):
        inputFile = args.input[i]
        sampleName = sampleNameList[i]
        binValDict = RunMetagene(inputFile, annoBedDict, args, kwargs)
        for binCoord in binValDict.keys():
            if binCoord not in binSampleValDict:
                binSampleValDict[binCoord] = defaultdict(dict)
            if 'feature' not in binSampleValDict[binCoord]:
                binSampleValDict[binCoord]['feature'] = binValDict[binCoord]['feature']
            binSampleValDict[binCoord][sampleName] = binValDict[binCoord]['value']
    ## output data
    with open(args.output, 'w') as out:
        ## print header
        row = ['feature', 'bin']
        row.extend(sampleNameList)
        out.write('\t'.join(row) + '\n')
        ## print data content
        for binCoord in sorted(binSampleValDict.keys()):
            feature = binSampleValDict[binCoord]['feature']
            row = [feature, str(binCoord)]
            valueRow = list(map(lambda x:str(binSampleValDict[binCoord][x]), sampleNameList))
            row.extend(valueRow)
            out.write('\t'.join(row) + '\n')
