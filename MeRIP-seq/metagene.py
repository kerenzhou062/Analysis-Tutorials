#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
import pybedtools
import subprocess
import tempfile
from multiprocessing import Pool
import gc
# in-house module
import bedutils

parser = argparse.ArgumentParser(
    description="This script is used for ploting metagene pattern from bed or bam based on gene annotations in bed12 format",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a', '--anno', action='store', type=str,
                    required=True,
                    help='The annotation bed12 file')
parser.add_argument('-b', '--bin', action='store', type=str,
                    choices=['average', 'constant'],
                    default='constant',
                    help='The type of bin size')
parser.add_argument('-c', '--cpu', action='store', type=int,
                    default=1,
                    help='Number of cpus to run this program (work for multiple input of --bed or --bam)')
parser.add_argument('-e', '--extend', action='store', type=int,
                    default=0,
                    help='Extend --extend bp around peak center')
parser.add_argument('-f', '--feature', nargs='+', type=str,
                    choices=["coding", "utr5", "cds", "utr3", "exon", "intron", "full"],
                    default='coding',
                    help='The bin features [coding:utr5,cds,utr3]')
parser.add_argument('-g', '--gene', nargs='+', type=str,
                    choices=["protein_coding", "non_coding", "all"],
                    default='protein_coding',
                    help='The bin features')
parser.add_argument('-i', '--bed', nargs='+', type=str,
                    help='The input bed files (bed3, bed6, bed12)')
parser.add_argument('-k', '--bam', nargs='+', type=str,
                    help='The input bam files (if --bed set, should be equal to --bed)')
parser.add_argument('-l', '--library', action='store', type=str,
                    choices=['unstranded', 'reverse', 'forward'],
                    default='reverse',
                    help='The library tye of bam files')
parser.add_argument('-m', '--method', action='store', type=str,
                    choices=['center', 'interval', 'exon'],
                    default='center',
                    help='The method for calculating bin coverage')
parser.add_argument('-n', '--name', nargs='+', type=str,
                    help='The sample name for --bed files')
parser.add_argument('-o', '--output', action='store', type=str,
                    default="bin.txt",
                    help='The output file')
parser.add_argument('-p', '--totalnum', action='store', type=int,
                    help='Set the total number of peaks for --bed or reads for --bam')
parser.add_argument('-r', '--memory', action='store', type=int,
                    default=5,
                    help='Set the memory (G) for sorting bed')
parser.add_argument('-s', '--smooth', action='store', type=str,
                    choices=['none', 'move', 'average'],
                    default='move',
                    help='The method for smoothing the curve')
parser.add_argument('-t', '--type', action='store', type=str,
                    choices=['density', 'percentage', 'number'],
                    default='density',
                    help='The type of output values (always "density" for --bam)')
parser.add_argument('-w', '--width', action='store', type=int,
                    default=5,
                    help='The span for smooth')
parser.add_argument('-x', '--temp', action='store', type=str,
                    default='/tmp',
                    help='The temporay directory')
parser.add_argument('-z', '--size', action='store', type=int,
                    default=100,
                    help='The bin size of each feature')
parser.add_argument('--matchid', action='store_true',
                    default=False,
                    help='Match input bed name in annotation bed12')
parser.add_argument('--deltmp', action='store_true',
                    default=False,
                    help='Delete ALL files matching "pybedtools.*.tmp" in the temp dir before running program')
parser.add_argument('--multiple', action='store_true',
                    default=False,
                    help='Use all multiple-aligned reads (use only one of them by default)')
parser.add_argument('--paired', action='store_true',
                    default=False,
                    help='Paired-end for bam files in --bam')
parser.add_argument('--reverse', action='store_true',
                    default=False,
                    help='Only report overlaps on the reverse strand (ignore by --library)')
parser.add_argument('--strand', action='store_true',
                    default=False,
                    help='Only report overlaps on the same strand (ignore by --library)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def FileExist(fList, argument):
    for file in fList:
        if os.path.isfile(file) is False:
            sys.exit('File "{0}" is not existed in --{1}!'.format(file, argument))

def ExonLenSum(blist):
    return sum(map(lambda x: x[1] - x[0], blist))

def SysSubCall(command):
    subprocess.call(command, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

def GetSampleName(fileLsit):
    sampleNameList = list()
    for file in fileLsit:
        basename = os.path.splitext(os.path.split(file)[-1])[0]
        sampleNameList.append(basename)
    return sampleNameList

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
                bedinfo.score = 1
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
                    bedLineRow.append(lineRow)
                lineNum += 1
    peakBed = pybedtools.BedTool(bedLineRow).sort()
    bedDict = {'bedtool':peakBed, 'totalNum':lineNum, 'source':'bed'}
    return bedDict

def BamToBed(bam, peakBed, destFile, args):
    bamToBedTemp = tempfile.NamedTemporaryFile(suffix='.tmp', prefix='pybedtools.tempfile', dir=args.temp, delete=True)
    if args.paired is True:
        fbamTemp = tempfile.NamedTemporaryFile(suffix='.tmp', prefix='pybedtools.tempfile', dir=args.temp, delete=True)
        command = 'samtools view -bf 65 {0} > {1}'.format(bam, fbamTemp.name)
        SysSubCall(command)
        command = 'bedtools bamtobed -i {0} > {1}'.format(fbamTemp.name, bamToBedTemp.name)
        SysSubCall(command)
        fbamTemp.close()
    else:
        command = 'bedtools bamtobed -i {0} > {1}'.format(bam, bamToBedTemp.name)
        SysSubCall(command)
    if args.multiple:
        unifyTempFile = tempfile.NamedTemporaryFile(suffix='.tmp', prefix='pybedtools.tempfile', dir=args.temp, delete=True)
        count = 1
        with open(bamToBedTemp.name, 'r') as f, open(unifyTempFile.name, 'w') as out:
            for line in f:
                row = line.split('\t')
                row[3] = '|'.join([row[3], str(count)])
                out.write('\t'.join(row))
                count += 1
        bamToBedTemp.close()
        command = 'sort -T {0} -S {1}G -k1,1 -k2,2n {2} -o {3}'.format(args.temp, args.memory, unifyTempFile.name, destFile)
        SysSubCall(command)
        unifyTempFile.close()
    else:
        command = 'sort -T {0} -S {1}G -k1,1 -k2,2n {2} -o {3}'.format(args.temp, args.memory, bamToBedTemp.name, destFile)
        SysSubCall(command)
        bamToBedTemp.close()
    ## get total read number
    command = 'wc -l {0}'.format(destFile)
    totalReadNum = int(bytes.decode(subprocess.check_output(command, shell=True)).split(' ')[0])
    ## construct bed
    bamToBed = pybedtools.BedTool(destFile)
    if peakBed is not None:
        kwargs = {'nonamecheck':True, 'u':True, 'sorted':True, 'S':True, 's':False}
        if args.library == 'unstranded':
            kwargs['S'] = False
        elif args.library == 'forward':
            kwargs['S'] = False
            kwargs['s'] = True
        bamToBed = bamToBed.intersect(peakBed, **kwargs).moveto(destFile)
        bamToBed = pybedtools.BedTool(destFile)
    bedDict = {'bedtool':bamToBed, 'totalNum':totalReadNum, 'source':'bam'}
    return bedDict

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
                bedLineRow.append(list(map(str, row)))
                fLenPreSum += coord[1] - coord[0]
    annoBed = pybedtools.BedTool(bedLineRow).sort()
    annoBedDict = {'bedtool':annoBed, 'bed12':bed12Dict, 'bed6':bed6Dict, 'binsize':binsizeDict}
    return annoBedDict

def RunMetagene(inputBedDict, annoBedDict, args, kwargs):
    inputBed = inputBedDict['bedtool']
    totalNum = inputBedDict['totalNum']
    bedSource = inputBedDict['source']
    annoBed = annoBedDict['bedtool']
    bed12Dict = annoBedDict['bed12']
    bed6Dict = annoBedDict['bed6']
    binsizeDict = annoBedDict['binsize']
    if bool(args.totalnum):
        totalNum = args.totalnum
    ## intersect rebuild inputBed and annoBed
    intersect = inputBed.intersect(annoBed, **kwargs)
    ## to reduce the memory usage
    inputBed = None
    interStatsDict = defaultdict(dict)
    for item in intersect:
        overlapRow = item.fields[0:6]
        annoRow = item.fields[6:]
        if bedSource == 'bam':
            peakName = overlapRow[3]
        else:
            peakName = overlapRow[3].split('##')[0]
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
        if peakName not in interStatsDict:
            interStatsDict[peakName] = defaultdict(int)
        interStatsDict[peakName][annoName] += interLen
    ## determin unique peakName-annoName relationship
    ## intersection lengh -> longest transcript
    sizeKey = args.feature + '_size'
    if args.feature == 'coding':
        sizeKey = 'exon_size'
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
                    preAnnoLen = bed12Dict[preAnnoName][sizeKey]
                    annoLen = bed12Dict[annoName][sizeKey]
                    if annoLen > preAnnoLen:
                        peakAnnoPairDict[peakName] = annoName
    del interStatsDict
    gc.collect()
    ## determin bin value
    binValDict = defaultdict(dict)
    for feature in binsizeDict.keys():
        binsize = binsizeDict[feature]
        binValDict[feature] = defaultdict(dict)
        for i in range(binsize):
            binValDict[feature][i]['sum'] = 0
            binValDict[feature][i]['peak'] = defaultdict(int)
    for item in intersect:
        overlapRow = item.fields[0:6]
        annoRow = item.fields[6:]
        if bedSource == 'bam':
            peakName = overlapRow[3]
        else:
            peakName = overlapRow[3].split('##')[0]
        uniqFeatureName = annoRow[3]
        annoName = annoRow[3].split('##')[0]
        if annoName == peakAnnoPairDict[peakName]:
            ## decode bed
            interStart = int(overlapRow[1])
            interEnd = int(overlapRow[2])
            interLen = interEnd - interStart
            feature, fLenPreSum = uniqFeatureName.split('##')[1:3]
            fLenPreSum = int(fLenPreSum)
            annoStart = bed6Dict[uniqFeatureName][1]
            annoEnd = bed6Dict[uniqFeatureName][2]
            annoStrand = bed6Dict[uniqFeatureName][5]
            mapLenPerbp = bed6Dict[uniqFeatureName][4]
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
    ## to reduce memory cost
    ## Deletes all temp files from the current session
    pybedtools.cleanup(verbose=False, remove_all=False)
    del intersect, peakAnnoPairDict
    gc.collect()
    ## generate final bin-value dict
    finalBinValDict = defaultdict(dict)
    featureList = [args.feature]
    count = 1
    if args.feature == 'coding':
        featureList = ['utr5', 'cds', 'utr3']
    for feature in featureList:
        binValList = list()
        for binCoord in sorted(binValDict[feature].keys()):
            binSum = binValDict[feature][binCoord]['sum']
            if bedSource == 'bam':
                binVal = binSum / totalNum
            else:
                if args.type == 'density':
                    binPeakNum = len(binValDict[feature][binCoord]['peak'].keys())
                    binVal = binPeakNum / totalNum * 100
                elif args.type == 'percentage':
                    binVal = binSum / totalNum * 100
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

def MultiThreadRun(index, iboolDict, annoBedDict, args, kwargs):
    if iboolDict['bam']:
        tempFile = tempfile.NamedTemporaryFile(suffix='.tmp', prefix='pybedtools.tempfile', dir=args.temp, delete=True)
    if iboolDict['both']:
        peakBed = pybedtools.BedTool(args.bed[index]).sort()
        bamFile = args.bam[index]
        inputBedDict = BamToBed(bamFile, peakBed, tempFile.name, args)
    elif iboolDict['bed']:
        inputBedDict = RebuildBed(args.bed[index], args.method, args.extend)
    else:
        peakBed = None
        bamFile = args.bam[index]
        inputBedDict = BamToBed(bamFile, peakBed, tempFile.name, args)
    ## retrieve bin-value relationships
    sampleName = args.name[index]
    binValDict = RunMetagene(inputBedDict, annoBedDict, args, kwargs)
    ## Deletes all temp files from the current session
    pybedtools.cleanup(verbose=False, remove_all=False)
    if iboolDict['bam']:
        tempFile.close()
    return [sampleName, binValDict]

# main program
if __name__ == '__main__':
    ## setting temporary dir for pybedtools
    pybedtools.set_tempdir(args.temp)
    if args.deltmp:
        pybedtools.cleanup(verbose=False, remove_all=True)
    ## judge arguments
    FileExist([args.anno], 'anno')
    if args.reverse and args.strand:
        sys.exit('--reverse and --strand should not be True at the same time!')
    if args.gene != 'protein_coding' and args.feature in ['coding', 'utr5', 'cds', 'utr3']:
        sys.exit('--feature should not be ["coding", "utr5", "cds", "utr3"] when --gene is not "protein_coding"!')
    ## bam and bed judgement
    iboolDict = defaultdict(bool)
    ibool = 0
    if bool(args.bed):
        FileExist(args.bed, 'bed')
        ibool += 1
        iboolDict['bed'] = True
    if bool(args.bam):
        FileExist(args.bam, 'bam')
        ibool += 1
        iboolDict['bam'] = True
    if ibool == 0:
        sys.exit('--bed or --bam should be set!')
    elif ibool == 2:
        iboolDict['both'] = True
    if iboolDict['both']:
        if len(args.bam) != len(args.bed):
            sys.exit('--bam and --bed should be matched!')
    ## get sample name list
    if bool(args.name):
        if iboolDict['bed']:
            if len(args.bed) != len(args.name):
                sys.exit('--bed and --name should be matched!')
        else:
            if len(args.bam) != len(args.name):
                sys.exit('--bam and --name should be matched!')
    else:
        if iboolDict['bed']:
            args.name = GetSampleName(args.bed)
        else:
            args.name = GetSampleName(args.bam)
    ## run programs
    kwargs = {'nonamecheck':True, 'wb':True, "sorted":True, 's':args.strand, 'S':args.reverse}
    if args.library == 'reverse' and iboolDict['bam']:
        kwargs['S'] = True
        kwargs['s'] = False
    elif args.library == 'unstranded' and iboolDict['bam']:
        kwargs['s'] = False
    ## construct annoBed from bed12
    annoBedDict = AnnoBed12ToBed6(args.anno, args.gene, args.feature, args.bin, args.size)
    ## multi-thread start
    pool = Pool(processes=args.cpu)
    resultList = []
    for i in range(len(args.name)):
        result = pool.apply_async(MultiThreadRun, args=(i, iboolDict, annoBedDict, args, kwargs,))
        resultList.append(result)
    pool.close()
    pool.join()
    ## multi-thread end
    binSampleValDict = defaultdict(dict)
    for result in resultList:
        sampleName, binValDict = result.get()
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
        row.extend(args.name)
        out.write('\t'.join(row) + '\n')
        ## print data content
        for binCoord in sorted(binSampleValDict.keys()):
            feature = binSampleValDict[binCoord]['feature']
            row = [feature, str(binCoord)]
            valueRow = list(map(lambda x:str(binSampleValDict[binCoord][x]), args.name))
            row.extend(valueRow)
            out.write('\t'.join(row) + '\n')
