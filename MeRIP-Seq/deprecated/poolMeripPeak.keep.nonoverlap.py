#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
import tempfile
import subprocess
from multiprocessing import Pool
from itertools import combinations
from copy import copy
import statistics
from glob import glob
import random
## custom modules
import bedutils

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='This script is used for pooling bed12 outputs from exomePeak|MeTPeak|MoAIMS, and re-calculate \
    the fold_enrichment and p-value. It only works with at least 2 biological replictes with bam.')
parser.add_argument('-a', '--anno', action='store', type=str,
                    help='Transcript annotations in bed12 format \
                    (4th column [gene_id:gene_name:gene_type:tx_id:tx_name:tx_type])). \
                    Required by --ip and --input.')
parser.add_argument('-e','--software', nargs='+', type=str,
                    required=True,
                    help='Software used for calling peaks for --peak (exomePeak|exomePeak2|MeTPeak|MACS2|SICER2|MoAIMS)')
parser.add_argument('-f', '--filter', action='store', type=float,
                    default=10,
                    help='IP reads number in peaks less than # will not be output to *.sigfc.bed')
parser.add_argument('-l', '--library', action='store', type=str,
                    choices=['unstranded', 'fr-firstrand','fr-secondstrand'],
                    default='fr-firstrand',
                    help='Library protocols of bam')
parser.add_argument('-n', '--name', nargs='+', type=str,
                    required=True,
                    help='Names for input --peak (eg. rep1 rep2)')
parser.add_argument('-o', '--output', action='store', type=str,
                    required=True,
                    help='The output directory')
parser.add_argument('-p', '--peak', nargs='+', type=str,
                    required=True,
                    help='M6A peak.bed (in bed12 format, 4thCol should be geneid) \
                    from (exomePeak|exomePeak2|MeTPeak|MACS2|SICER2|MoAIMS, refineBed12.py)')
parser.add_argument('-r', '--random', action='store', type=int,
                    default=100,
                    help='randomly pool the final peaks to avoid any overlaps')
parser.add_argument('-s', '--shrink', action='store', type=str,
                    choices=['none', 'apeglm', 'ashr'],
                    default='none',
                    help='The method used for shrinking the fold change (lfcShrink() in DESeq2)')
parser.add_argument('-t', '--test', action='store', type=str,
                    choices=['Wald', 'LRT'],
                    default='Wald',
                    help='The test method for calculating p-value (DESeq2)')
parser.add_argument('--extsize', action='store', type=int,
                    default=0,
                    help='extend reads in 5\'->3\' direction to # bp (only works in single-end data)')
parser.add_argument('--fraction', action='store', type=float,
                    default=0.1,
                    help='The overlap fraction of peaks should be more larger # when overlapping peaks')
parser.add_argument('--ip', nargs='*', type=str,
                    help='Alignment bam files of conrtol IP samples')
parser.add_argument('--input', nargs='*', type=str,
                    help='Alignment bam files of conrtol input samples')
parser.add_argument('--fragmentSize', action='store', type=int,
                    default=600,
                    help='The fragment size of library, important for generating non-peak region')
parser.add_argument('--log2fc', action='store', type=float,
                    default=1,
                    help='Cutoff log2FoldChange for significant peaks (only used when --ip and --input are set, DESeq2)')
parser.add_argument('--max', action='store', type=int,
                    default=5000000,
                    help='Filter peaks larger than --max')
parser.add_argument('--min', action='store', type=int,
                    default=20,
                    help='Filter peaks smaller than --min (the overlap peaks should be larger than #bp when overlapping peaks)')
parser.add_argument('--prefix', action='store', type=str,
                    default='m6a_peak',
                    help='File prefix for output results')
parser.add_argument('--pval', action='store', type=float,
                    default=0.05,
                    help='Cutoff of pval for significant peaks (only used when --ip and --input are set, DESeq2)')
parser.add_argument('--padj', action='store', type=float,
                    default=0.1,
                    help='Cutoff of pvalAdjust for significant peaks (only used when --ip and --input are set, DESeq2)')
parser.add_argument('--padjMethod', action='store', type=str,
                    choices=["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"],
                    default='BH',
                    help='The method to calculate the padj (in R, ?p.adjust)')
parser.add_argument('--thread', action='store', type=int,
                    default=1,
                    help='The number of threads used to run in parallel (will greatly reduce running time)')
parser.add_argument('--tempDir', action='store', type=str,
                    help='directory for storing temporary files, default is --output')
parser.add_argument('--nofraction', action='store_true',
                    default=False,
                    help='Do not count multiple alignment or overlapped reads as fraction (disable "--fraction" in featureCounts)')
parser.add_argument('--pairend', action='store_true',
                    default=False,
                    help='The input bams are pair-end')
parser.add_argument('--keepTemp', action='store_true',
                    default=False,
                    help='Keep the temporary files')
parser.add_argument('--unique', action='store_true',
                    default=False,
                    help='Counts uniquely aligned and overlapped reads or fragments only')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def DeleteFile(file):
    try:
        os.remove(file)
    except:
        pass

def GetSortPeak(peakRowList):
    peakIndexDict = {}
    for i in range(len(peakRowList)):
        peakRowList[i][3] = '=='.join([peakRowList[i][3], str(i)])
        peakIndexDict[peakRowList[i][3]] = i
    sortedPeakRowList = sorted(peakRowList, key=lambda x:(x[0], int(x[1]), int(x[2])))
    ## get sorted index map to original index
    indexMapDict = {}
    for i in range(len(sortedPeakRowList)):
        originalIndex = peakIndexDict[sortedPeakRowList[i][3]]
        indexMapDict[i] = originalIndex
        sortedPeakRowList[i][3] = sortedPeakRowList[i][3].split('==')[0]
    return [sortedPeakRowList, indexMapDict]

def SubprocessToList(command, errorFlag=False):
    if errorFlag is True:
        runResList = bytes.decode(subprocess.check_output(command, shell=True, stderr=None)).split('\n')
    else:
        runResList = bytes.decode(subprocess.check_output(command, shell=True, stderr=subprocess.PIPE)).split('\n')
    runResRowList = list(map(lambda x:x.split('\t'), list(filter(lambda x: bool(x), runResList))))
    return runResRowList

def GetReadLength(bamFile):
    command = 'samtools view {} | head -n 20 | awk \'{{print $10}}\''.format(bamFile)
    readRowList = SubprocessToList(command)
    readLenList = []
    for readRow in readRowList:
        read = readRow[0]
        readLenList.append(len(read))
    readLength = int(statistics.mean(readLenList))
    return readLength

def FindNonPeakRegion(peakRowList, annoBed12, fragmentSize, tempDir):
    ## merge exons in annotations
    fragmentSize = fragmentSize * 2
    annoMergeExonTmp = tempfile.NamedTemporaryFile(suffix='.tmp', dir=tempDir, delete=True)
    command = "bedtools bed12tobed6 -i {} | sort -t$'\\t' -k1,1 -k2,2n -T {} | "
    command += "bedtools merge -i stdin -s -c 6 -o distinct | "
    command += 'awk \'BEGIN{{FS="\t";OFS="\t"}}{{print $1,$2,$3,"region:"FNR, 0, $4}}\' > {}'
    command = command.format(annoBed12, tempDir, annoMergeExonTmp.name)
    __ = SubprocessToList(command)
    ## construct peak file
    peakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', dir=tempDir, delete=True)
    peakMergeExonTmp = tempfile.NamedTemporaryFile(suffix='.tmp', dir=tempDir, delete=True)
    with open(peakTmp.name, 'w') as temp:
        for row in peakRowList:
            temp.write('\t'.join(row) + '\n')
    ## merge exons in peak
    command = "bedtools bed12tobed6 -i {} | sort -t$'\\t' -k1,1 -k2,2n -T {} | "
    command += "bedtools merge -i stdin -s -c 6 -o distinct | "
    command += 'awk \'BEGIN{{FS="\t";OFS="\t"}}{{print $1,$2,$3,"peak:"FNR, 0, $4}}\' > {}'
    command = command.format(peakTmp.name, tempDir, peakMergeExonTmp.name)
    __ = SubprocessToList(command)
    ## substract annoMergeExon by peakMergeExon, then overlap with peakMergeExonTmp
    ## distance between 2 features should be in fragmentSize long as least
    command = 'bedtools subtract -a {} -b {} -s | '
    command += 'awk -v size="{}" \'BEGIN{{FS="\t";OFS="\t"}}{{$2=$2-size;$3=$3+size; if($2<0){{$2=0}} print $1,$2,$3,"region:"FNR, 0, $6}}\' |'
    command += 'bedtools intersect -a stdin -b {} -s -v |'
    command += 'awk -v size="{}" \'BEGIN{{FS="\t";OFS="\t"}}{{$2=$2+size;$3=$3-size; if($3 - $2>0){{print $1,$2,$3,"region:"FNR, 0, $6}}}}\''
    command = command.format(annoMergeExonTmp.name, peakMergeExonTmp.name, fragmentSize, peakMergeExonTmp.name, fragmentSize)
    nonOverlapRegionRowList = SubprocessToList(command)
    ## substract annoMergeExon by peakMergeExon, then overlap with peakMergeExonTmp
    ## to get non-overlap regions that adjacent to the peaks
    command = 'bedtools subtract -a {} -b {} -s | '
    command += 'awk -v size="{}" \'BEGIN{{FS="\t";OFS="\t"}}{{$2=$2-size;$3=$3+size; if($2<0){{$2=0}} print $1,$2,$3,"region:"FNR, 0, $6}}\' |'
    command += 'bedtools intersect -a stdin -b {} -s -wa -wb'
    command = command.format(annoMergeExonTmp.name, peakMergeExonTmp.name, fragmentSize, peakMergeExonTmp.name)
    adjacentRegionRowList = SubprocessToList(command)
    ## delete temp file
    peakTmp.close()
    annoMergeExonTmp.close()
    peakMergeExonTmp.close()
    ## generate final regions that do no overlap with peaks
    ## to avoid regions overlap with reads that overlaping with peaks
    ## substract with the fragment size
    adjacentRegionDict = defaultdict(dict)
    for adjRegionRow in adjacentRegionRowList:
        annoRegionBed = bedutils.buildbed(adjRegionRow[0:6])
        peakRegionBed = bedutils.buildbed(adjRegionRow[6:])
        adjacentId = 'NonPeak:' + annoRegionBed.name
        if peakRegionBed.end - annoRegionBed.start <= fragmentSize:
            adjRegionRow[1] = annoRegionBed.start + fragmentSize * 2
            if adjacentId not in adjacentRegionDict:
                adjacentRegionDict[adjacentId] = adjRegionRow[0:6]
            else:
                adjacentRegionDict[adjacentId][1] = adjRegionRow[1]
        elif annoRegionBed.end - peakRegionBed.start <= fragmentSize:
            adjRegionRow[2] = annoRegionBed.end - fragmentSize * 2
            if adjacentId not in adjacentRegionDict:
                adjacentRegionDict[adjacentId] = adjRegionRow[0:6]
            else:
                adjacentRegionDict[adjacentId][2] = adjRegionRow[2]
    ## final adjacentRowList
    ## chr, start, end, name, score, strand
    adjacentRegionRowList = list(map(lambda x:adjacentRegionDict[x], sorted(adjacentRegionDict.keys())))
    ## final non-peak region
    nonPeakRegionRowList = nonOverlapRegionRowList + adjacentRegionRowList
    return nonPeakRegionRowList

def RegionToSafMatrix(peakRowList, nonPeakRegionRowList, safMatrixFile):
    safRowList = [['GeneID', 'Chr', 'Start', 'End', 'Strand']]
    ##  convert peakRow into the SAF-format
    for peakRow in peakRowList:
        peakBed = bedutils.buildbed(peakRow).decode()
        peakid = peakBed.name
        chrom = peakBed.chr
        strand = peakBed.strand
        for exon in peakBed.exon:
            safRow = [peakid, chrom, exon[0] + 1, exon[1], strand]
            safRowList.append(safRow)
    ## get exons without peak intersection
    ## this step aims to balance the final reads library to meet the hypothesis of DESeq2
    ## that is the most of genes are not changed among libraries
    count = 1
    for regionRow in nonPeakRegionRowList:
        if (int(regionRow[2]) - int(regionRow[1])) < args.min:
            continue
        uniqid = 'NonPeak:region:' + str(count)
        safRow = [uniqid, regionRow[0], int(regionRow[1]) + 1, regionRow[2], regionRow[5]]
        safRowList.append(safRow)
        count += 1
    ## construct final saf file library
    with open(safMatrixFile, 'w') as temp:
        for row in safRowList:
            temp.write('\t'.join(map(str, row)) + '\n')
    return safMatrixFile

def CollapseDupBedRow(bedPeakRowList, bedLabelRowList, peakNameSep):
    ## keep duplicated genomic coordinates (bed12) in one record between bedPeakRowList
    ## bedPeakRowList = [peakRowList1, peakRowList2, ...]
    ## bedLabelRowList = [labelRowList1, labelRowList2, ...]
    peakKeyDict = defaultdict(list)
    indexToKeyDict = defaultdict(dict)
    for i in range(len(bedPeakRowList)):
        peakRowList = bedPeakRowList[i]
        for j in range(len(peakRowList)):
            peakRow = copy(peakRowList[j])
            geneid = peakRow[3].split('|')[0].split(peakNameSep)[0]
            peakRow[3] = geneid
            peakRow[4] = '1'
            peakRow[6] = peakRow[6]
            peakRow[7] = peakRow[7]
            peakRow[8] = '0'
            peakRow[10] = peakRow[10].strip(',')
            peakRow[11] = peakRow[11].strip(',')
            key = '='.join(peakRow[0:12])
            ## record with key and index
            peakKeyDict[key].append([i,j])
            indexToKeyDict[i][j] = key
    ## construct new bedPeakRowList and new bedLabelRowList
    newBedPeakRowList = [[] for i in bedPeakRowList]
    newBedLabelRowList = [[] for i in bedLabelRowList]
    keyUsedDict = defaultdict(dict)
    for i in range(len(bedPeakRowList)):
        peakRowList = bedPeakRowList[i]
        labelRowList = bedLabelRowList[i]
        for j in range(len(peakRowList)):
            key = indexToKeyDict[i][j]
            if len(peakKeyDict[key]) == 1:
                newBedPeakRowList[i].append(peakRowList[j])
                newBedLabelRowList[i].append(labelRowList[j])
                keyUsedDict[key] = 1
            else:
                if key not in keyUsedDict:
                    keyUsedDict[key] = 1
                    nameList = []
                    softwareList = []
                    for indexGroup in peakKeyDict[key]:
                        labelList = bedLabelRowList[indexGroup[0]][indexGroup[1]]
                        nameList.extend(labelList[0].split(','))
                        softwareList.extend(labelList[1].split(','))
                    name = ','.join(sorted(set(nameList)))
                    software = ','.join(sorted(set(softwareList)))
                    labelRow = [name, software]
                    ## store to newBedPeakRowList and newBedLabelRowList
                    newBedPeakRowList[i].append(peakRowList[j])
                    newBedLabelRowList[i].append(labelRow)
    return [newBedPeakRowList, newBedLabelRowList]

def RunFeatureCounts(safMatrix, countFile, bamFiles, thread, libraryStrand, pairend, extsize, nofractionFlag, uniqueFlag, errorFlag=False):
    if uniqueFlag is True:
        uniqueParam = ""
    else:
        if nofractionFlag is True:
            uniqueParam = "-C -M -O"
        else:
            uniqueParam = "-C -M -O --fraction"
    if pairend is False:
        command = 'featureCounts -a {} -F SAF -o {} -s {} -T {} --readExtension3 {} {} {}'
        command = command.format(safMatrix, countFile, libraryStrand, thread, extsize, uniqueParam, bamFiles)
    else:
        command = 'featureCounts -a {} -F SAF -o {} -s {} -T {} -p --countReadPairs {} {}'
        command = command.format(safMatrix, countFile, libraryStrand, thread, uniqueParam, bamFiles)
    __ = SubprocessToList(command, errorFlag)
    ## get reads
    peakReadDict = defaultdict(dict)
    with open(countFile, 'r') as f:
         __ = f.readline()
         ## get header row and mark index of bamfile
         headerRow = f.readline().strip().split('\t')
         rowLength = len(headerRow)
         toIndexDict = defaultdict(dict)
         for i in range(6,rowLength):
             bamFile = headerRow[i]
             uniqBamIndex = bamDict['bamToIndex'][bamFile]
             toIndexDict[i] = uniqBamIndex
         ##  get reads in libraries
         for line in f:
             row = line.strip().split('\t')
             peakid = row[0]
             for i in range(6,rowLength):
                 uniqBamIndex = toIndexDict[i]
                 readsNum = row[i]
                 peakReadDict[peakid][uniqBamIndex] = readsNum
    ## delete file from featureCounts
    DeleteFile(countFile)
    DeleteFile(countFile + '.summary')
    ## return final reads dict
    return peakReadDict

def SeparatePeak(peakRowList, labelRowList):
    ## seperate input peakRowList into 2 parts:
    ## uniq (do not overlap with any peaks) and overlap part
    global peakNameSep
    peakNum = len(peakRowList)
    indexList = list(range(peakNum))
    genePeakDict = defaultdict(list)
    indexLabelDict = defaultdict(int)
    for i in indexList:
        geneid = peakRowList[i][3].split(peakNameSep)[0]
        genePeakDict[geneid].append(i)
        indexLabelDict[i] = 0
    ## label peak as uniq (==0) or overlap (>0)
    for geneid in sorted(genePeakDict.keys()):
        geneIndexList = genePeakDict[geneid]
        geneIndexNum = len(geneIndexList)
        if len(geneIndexList) > 1:
            for i in range(geneIndexNum):
                indexI = geneIndexList[i]
                iBedOps = bedutils.bed12ops(peakRowList[indexI])
                for j in range(i + 1, geneIndexNum):
                    indexJ = geneIndexList[j]
                    jBedOps = bedutils.bed12ops(peakRowList[indexJ])
                    if (iBedOps.a.bcount == 1 or jBedOps.a.bcount == 1):
                        intersect = iBedOps.intersect(peakRowList[indexJ], score='max', s=True, tx=True, part=True, cds=False)
                    else:
                        intersect = iBedOps.intersect(peakRowList[indexJ], score='max', s=True, tx=True, part=False, cds=False)
                    if bool(intersect) is True:
                        indexLabelDict[indexI] += 1
                        indexLabelDict[indexJ] += 1
    uniqIndexList = list(filter(lambda x: indexLabelDict[x] == 0, indexList))
    overlapIndexList = list(filter(lambda x: indexLabelDict[x] > 0, indexList))
    ## get separated list
    uniqPeakRowList = [peakRowList[i] for i in uniqIndexList]
    uniqLabelRowList = [labelRowList[i] for i in uniqIndexList]
    overlapPeakRowList = [peakRowList[i] for i in overlapIndexList]
    overlapLabelRowList = [labelRowList[i] for i in overlapIndexList]
    ## return results
    return [uniqPeakRowList, uniqLabelRowList, overlapPeakRowList, overlapLabelRowList]


def RandomPool(keyword, peakRowList, labelRowList, limit, finalFlag=True):
    global peakNameSep
    global softwareDict
    ## randomly split peakRowList into 2 equal list, and then perform PoolPeak
    ## to collapse overlapped peaks in the finalPoolPeakRowList
    interation = 0
    prevOverlapNum = 0
    sameFlag = 0
    randomStopCutoff = 5
    ## seperate peaks firstly, greatly speed up
    uniqPeakRowList, uniqLabelRowList, poolPeakRowList, poolLabelRowList = SeparatePeak(peakRowList, labelRowList)
    ## only run ramdom pooling on overlap peaks
    for seed in range(1, limit + 1):
        peakNum = len(poolPeakRowList)
        halfNum = int(peakNum / 2)
        indexList = list(range(peakNum))
        random.seed(seed)
        random.shuffle(indexList)
        index1List = indexList[0:halfNum]
        index2List = indexList[halfNum:]
        peak1RowList = [poolPeakRowList[i] for i in index1List]
        label1RowList = [poolLabelRowList[i] for i in index1List]
        peak2RowList = [poolPeakRowList[i] for i in index2List]
        label2RowList = [poolLabelRowList[i] for i in index2List]
        poolPeakRowList, poolLabelRowList = PoolPeak(peak1RowList, peak2RowList, label1RowList, label2RowList, finalFlag)
        __, __, testOverlapPeakRowList, __ = SeparatePeak(poolPeakRowList, poolLabelRowList)
        if len(testOverlapPeakRowList) == 0:
            interation += 1
            break
        else:
            ## if the overlap number always the same, then stop ramdon pool
            if sameFlag >= randomStopCutoff:
                interation += 1
                break
            if prevOverlapNum == len(testOverlapPeakRowList):
                sameFlag += 1
            else:
                prevOverlapNum = len(testOverlapPeakRowList)
                sameFlag = 0
        sys.stderr.write("{}: ramdom pool interation {}, {} overlaps remain\n".format(keyword, interation, len(testOverlapPeakRowList)))
        ## count total interation times
        interation += 1
    if finalFlag is True:
        sys.stderr.write("{}: ramdom pool finally stop at interation {} \n".format(keyword, interation))
    ## return results
    peakRowList = uniqPeakRowList + poolPeakRowList
    labelRowList = uniqLabelRowList + poolLabelRowList
    return [peakRowList, labelRowList]

def ExonToBed12Row(chrom, strand, name, score, rgb, exonList):
    start = min(map(lambda x:x[0], exonList))
    end = max(map(lambda x:x[1], exonList))
    blockSizesList = []
    blockStartsList = []
    for exon in exonList:
        blockSize = exon[1] - exon[0]
        blockStart = exon[0] - start
        blockSizesList.append(blockSize)
        blockStartsList.append(blockStart)
    ## joint the block size and block starts with ','
    blockSize = ','.join(map(str, blockSizesList))
    blockStarts = ','.join(map(str, blockStartsList))
    row = [chrom, start, end, name, score, strand, start, end, rgb, len(exonList), blockSize, blockStarts]
    row = list(map(str, row))
    return row

def GetNonOverlapBed(tbuildBed, cbuildBed, minLength, overhangFlag=True):
    ## return the remaining non-intersected exons of builBedA from cbuildBed in bed12 format
    ## intersected exon located at first or the last exon of A
    ## eg. A:a-b-c, B:c, then return a-b-c(non-intersected part)
    if cbuildBed.a.bcount == 1 and tbuildBed.a.bcount > 1:
        refExon = cbuildBed.a.exon[0]
        count = 0
        indexTag = -1
        for i in range(len(tbuildBed.a.exon)):
            exon = tbuildBed.a.exon[i]
            distance = min(refExon[1], exon[1]) - max(refExon[0], exon[0])
            if distance > 0:
                count += 1
                indexTag = i
        if count == 1 and (indexTag == 0 or indexTag == len(tbuildBed.a.exon) - 1):
            remainExonList = []
            chrom = tbuildBed.a.chr
            strand = tbuildBed.a.strand
            name = tbuildBed.a.name
            score = tbuildBed.a.score
            rgb = tbuildBed.a.rgb
            for i in range(len(tbuildBed.a.exon)):
                exon = tbuildBed.a.exon[i]
                if i == indexTag:
                    refBedRow = ['chr1', refExon[0], refExon[1], 'name1', 0, '+']
                    exonBedRow = ['chr1', exon[0], exon[1], 'name2', 0, '+']
                    intersect = bedutils.bed6ops(refBedRow).intersect(exonBedRow)
                    if overhangFlag is False:
                        ## overhang-overlap is not allowed for returnning non-intersected part
                        if refExon[0] < exon[0] or refExon[1] > exon[1]:
                            return False
                    if i == 0:
                        remainExon = [intersect.a.end + minLength, exon[1]]
                    else:
                        remainExon = [exon[0], intersect.a.start - minLength]
                    if remainExon[1] - remainExon[0] >= minLength:
                        remainExonList.append(remainExon)
                else:
                    remainExonList.append(exon)
            bedRow = ExonToBed12Row(chrom, strand, name, score, rgb, remainExonList)
            buildbed = bedutils.buildbed(bedRow)
            return buildbed
    else:
        return False

def IntersectPeakPool(cpeakRow, tpeakRow, peakLabelDict, minLength, recordPeakDict, remainPeakDict, finalFlag=False, remainFlag=False):
    ## remainFlag, the cpeakRow and tpeakRow 
    ## intersect 2 peaks
    global fractionCutoff
    nonOverlapPerCutoff = 0.7
    cpeakid = cpeakRow[3]
    tpeakid = tpeakRow[3]
    ## create bedutils object
    cpeakBedOps = bedutils.bed12ops(cpeakRow)
    tpeakBedOps = bedutils.bed12ops(tpeakRow)
    ## determine the software priority
    cpeakSoftwareList = peakLabelDict[cpeakid][1].split(',')
    cpeakSoftPriv = min(map(lambda x:softwareDict[x], cpeakSoftwareList))
    tpeakSoftwareList = peakLabelDict[tpeakid][1].split(',')
    tpeakSoftPriv = min(map(lambda x:softwareDict[x], tpeakSoftwareList))
    ## start to intersect
    if cpeakBedOps.a.bcount == 1 or tpeakBedOps.a.bcount == 1:
        intersect = cpeakBedOps.intersect(tpeakRow, score='max', s=True, tx=True, part=True, cds=False)
    else:
        intersect = cpeakBedOps.intersect(tpeakRow, score='max', s=True, tx=True, part=False, cds=False)
    ## if intersect
    if bool(intersect) is True:
        overlapFractionBool = ( intersect.a.exonlength/cpeakBedOps.a.exonlength >= fractionCutoff or intersect.a.exonlength/tpeakBedOps.a.exonlength )
        if intersect.a.exonlength >= minLength and overlapFractionBool is True:
            overlapBool = True
        else:
            overlapBool = False
    else:
        overlapBool = False
    if overlapBool is True:
        if remainFlag is False:
            if 'valid' not in recordPeakDict[cpeakid]:
                recordPeakDict[cpeakid]['valid'] = 1
            else:
                recordPeakDict[cpeakid]['valid'] += 1
            if 'valid' not in recordPeakDict[tpeakid]:
                recordPeakDict[tpeakid]['valid'] = 1
            else:
                recordPeakDict[tpeakid]['valid'] += 1
        bcountBool = (cpeakBedOps.a.bcount == 1 and tpeakBedOps.a.bcount > 1) or (tpeakBedOps.a.bcount == 1 and cpeakBedOps.a.bcount > 1)
        if cpeakSoftPriv == tpeakSoftPriv:
            poolPeakBedOps = intersect
            if bcountBool:
                ## return non-intersected region if possible
                if cpeakBedOps.a.bcount > 1:
                    largerBedOps, smallBedOps = cpeakBedOps, tpeakBedOps
                else:
                    largerBedOps, smallBedOps = tpeakBedOps, cpeakBedOps
                poolPeakBedOps = largerBedOps
                if finalFlag is False and remainFlag is False:
                    ## get non overlap part of tpeak
                    peakNonOverlapBed = GetNonOverlapBed(largerBedOps, smallBedOps, minLength, False)
                    if bool(peakNonOverlapBed) is True:
                        nonOverlapPer = peakNonOverlapBed.exonlength / (intersect.a.exonlength + minLength + peakNonOverlapBed.exonlength)
                        if nonOverlapPer >= nonOverlapPerCutoff:
                            remainPeakDict[largerBedOps.a.name].append(peakNonOverlapBed.list)
                    else:
                        poolPeakBedOps = largerBedOps
                        poolPeakBedOps.a.list[3] = intersect.a.name
        elif cpeakSoftPriv < tpeakSoftPriv:
            poolPeakBedOps = cpeakBedOps
            poolPeakBedOps.a.list[3] = intersect.a.name
            if finalFlag is False and remainFlag is False:
                ## get non overlap part of tpeak if cpeak contains only 1 exon and tpeak contains 2 or more exons
                peakNonOverlapBed = GetNonOverlapBed(tpeakBedOps, cpeakBedOps, minLength, True)
                if bool(peakNonOverlapBed) is True:
                    nonOverlapPer = peakNonOverlapBed.exonlength / (intersect.a.exonlength + minLength + peakNonOverlapBed.exonlength)
                    if nonOverlapPer >= nonOverlapPerCutoff:
                        remainPeakDict[tpeakid].append(peakNonOverlapBed.list)
        else:
            poolPeakBedOps = tpeakBedOps
            poolPeakBedOps.a.list[3] = intersect.a.name
            if finalFlag is False and remainFlag is False:
                ## get non overlap part of cpeak if tpeak contains only 1 exon and cpeak contains 2 or more exons
                peakNonOverlapBed = GetNonOverlapBed(cpeakBedOps, tpeakBedOps, minLength, True)
                if bool(peakNonOverlapBed) is True:
                    nonOverlapPer = peakNonOverlapBed.exonlength / (intersect.a.exonlength + minLength + peakNonOverlapBed.exonlength)
                    if nonOverlapPer >= nonOverlapPerCutoff:
                        remainPeakDict[cpeakid].append(peakNonOverlapBed.list)
        ## record the intersected peak
        poolPeakRow = poolPeakBedOps.a.list
    else:
        poolPeakRow = []
        if remainFlag is False:
            recordPeakDict[cpeakid]['invalid'] = 1
            recordPeakDict[tpeakid]['invalid'] = 1
    ## final results
    if remainFlag is True:
        return poolPeakRow
    else:
        return [poolPeakRow, recordPeakDict, remainPeakDict]

def PoolPeak(cpeakRowList, tpeakRowList, clabelRowList, tlabelRowList, finalFlag=False):
    ## label1RowList: [ ['shCont_rep1', 'exomePeak'], ['shCont_rep1', 'exomePeak'] ...]
    ## make peak id unique
    global peakNameSep
    global softwareDict
    global minLength
    global maxLength
    global tempDir
    cpeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', dir=tempDir, delete=True)
    peakRowDict = defaultdict(dict)
    peakLabelDict = defaultdict(dict)
    count = 1
    with open(cpeakTmp.name, 'w') as out:
        for i in range(len(cpeakRowList)):
            row = cpeakRowList[i]
            row[3] = peakNameSep.join([row[3], str(count)])
            ## record peakid -> row
            peakRowDict[row[3]] = row
            ## record peakid -> labelRow ([name, software])
            labelRow = clabelRowList[i]
            peakLabelDict[row[3]] = labelRow
            ## write to temp file
            out.write('\t'.join(row) + '\n')
            count += 1
    tpeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', dir=tempDir, delete=True)
    with open(tpeakTmp.name, 'w') as out:
        for i in range(len(tpeakRowList)):
            row = tpeakRowList[i]
            row[3] = peakNameSep.join([row[3], str(count)])
            ## record peakid -> row
            peakRowDict[row[3]] = row
            ## record peakid -> labelRow
            labelRow = tlabelRowList[i]
            peakLabelDict[row[3]] = labelRow
            ## write to temp file
            out.write('\t'.join(row) + '\n')
            count += 1
    ## get peaks only found in control groups
    command = 'bedtools intersect -a {} -b {} -v -s -split -nonamecheck'.format(cpeakTmp.name, tpeakTmp.name)
    cuniqPeakRowList = SubprocessToList(command)
    ## get peaks only found in treated groups
    command = 'bedtools intersect -a {} -b {} -v -s -split -nonamecheck'.format(tpeakTmp.name, cpeakTmp.name)
    tuniqPeakRowList = SubprocessToList(command)
    ## get overlap peaks in control and treated groups
    command = 'bedtools intersect -a {} -b {} -split -s -wa -wb -nonamecheck'.format(cpeakTmp.name, tpeakTmp.name)
    overlapPeakRowList = SubprocessToList(command)
    ## close temp file
    cpeakTmp.close()
    tpeakTmp.close()
    ## check whether there are any overlap peaks
    poolPeakRowList = []
    recordPeakDict = defaultdict(dict)
    remainPeakDict = defaultdict(list)
    if len(overlapPeakRowList) > 0:
        intersectDict = defaultdict(list)
        for row in overlapPeakRowList:
            cpeakid = row[3]
            tpeakid = row[15]
            cpeakGeneName = cpeakid.split(peakNameSep)[0]
            tpeakGeneName = tpeakid.split(peakNameSep)[0]
            ## filter out the overlapped peaks on different genes
            if cpeakGeneName != tpeakGeneName:
                recordPeakDict[cpeakid]['invalid'] = 1
                recordPeakDict[tpeakid]['invalid'] = 1
            else:
                intersectDict[cpeakid].append(tpeakid)
        ## start to pool peaks
        for cpeakid in sorted(intersectDict.keys()):
            cpeakRow = peakRowDict[cpeakid]
            tpeakidList = intersectDict[cpeakid]
            for tpeakid in tpeakidList:
                tpeakRow = peakRowDict[tpeakid]
                poolPeakRow, recordPeakDict, remainPeakDict = IntersectPeakPool(cpeakRow, tpeakRow, peakLabelDict, minLength, recordPeakDict, remainPeakDict, finalFlag, False)
                if bool(poolPeakRow):
                    poolPeakRowList.append(poolPeakRow)
    peakRowNonOverlapList = []
    ## deal with remain peak-part
    for peakid in sorted(remainPeakDict.keys()):
        recordFlagDict = defaultdict(dict)
        remainPeakRowList, __ = GetSortPeak(remainPeakDict[peakid])
        remainNum = len(remainPeakRowList)
        if remainNum == 1:
            peakRowNonOverlapList.append(remainPeakRowList[0])
        else:
            singleBlockIndexList = []
            multiBlockIndexList = []
            for i in range(remainNum):
                cpeakBed = bedutils.buildbed(remainPeakRowList[i])
                if cpeakBed.bcount > 1:
                    multiBlockIndexList.append(i)
                else:
                    singleBlockIndexList.append(i)
            for indexBlockList in [multiBlockIndexList, singleBlockIndexList]:
                if len(indexBlockList) == 1:
                    indexI = indexBlockList[0]
                    peakRowNonOverlapList.append(remainPeakRowList[indexI])
                    continue
                for i in range(len(indexBlockList)):
                    indexI = indexBlockList[i]
                    cpeakRow = remainPeakRowList[indexI]
                    if indexI in recordFlagDict:
                        continue
                    for j in range(i+1, len(indexBlockList)):
                        indexJ = indexBlockList[j]
                        if indexJ in recordFlagDict:
                            continue
                        poolPeakRow = IntersectPeakPool(cpeakRow, remainPeakRowList[indexJ], peakLabelDict, minLength, 'temp', 'temp', True, True)
                        if bool(poolPeakRow):
                            poolPeakRow[3] = poolPeakRow[3].split('|')[0]
                            cpeakRow = poolPeakRow
                            recordFlagDict[i] = 1
                            recordFlagDict[j] = 1
                    peakRowNonOverlapList.append(cpeakRow)
    # retrive non intersect peaks
    for peakid in sorted(recordPeakDict.keys()):
        if 'invalid' in recordPeakDict[peakid]:
            if 'valid' not in recordPeakDict[peakid]:
                peakRowNonOverlapList.append(peakRowDict[peakid])
            else:
                if recordPeakDict[peakid]['valid'] < 1:
                    peakRowNonOverlapList.append(peakRowDict[peakid])
    ## pool the unique and merged bed rows
    peakRowList = cuniqPeakRowList + tuniqPeakRowList + peakRowNonOverlapList + poolPeakRowList
    ## sort rows by chromosomes and start positions
    peakRowList, __ = GetSortPeak(peakRowList)
    ### collapse duplicated records based on genomic positions
    labelRowList = []
    for i in range(len(peakRowList)):
        peakRow = copy(peakRowList[i])
        peakNameRow = peakRow[3].split('|')
        ## get label name and software, [[rep1, rep2], ['rep1']] -> ['rep1', 'rep2']
        nameList = list(set([x for y in map(lambda x:peakLabelDict[x][0].split(','), peakNameRow) for x in y]))
        softwareList = list(set([x for y in map(lambda x:peakLabelDict[x][1].split(','), peakNameRow) for x in y]))
        labelRow = [','.join(sorted(nameList)), ','.join(sorted(softwareList))]
        labelRowList.append(labelRow)
    [peakRowList], [labelRowList] = CollapseDupBedRow([peakRowList], [labelRowList], peakNameSep)
    ## construct final peak row list and label row list
    finalPeakRowList = []
    finalLabelRowList = []
    for i in range(len(peakRowList)):
        peakRow = copy(peakRowList[i])
        peakRow[3] = peakRow[3].split(peakNameSep)[0]
        peakBed = bedutils.buildbed(peakRow)
        labelRow = labelRowList[i]
        ## filter out peaks that larger or smaller than required
        if peakBed.exonlength < minLength or peakBed.exonlength > maxLength:
            continue
        finalPeakRowList.append(peakRow)
        finalLabelRowList.append(labelRow)
    return [finalPeakRowList, finalLabelRowList]

def RunDESeq2(poolPeakRowList, fragmentSize, bamFileList, libraryStrand, finalFlag=False):
    global bamDict
    global geneAnnoFile
    global tempDir
    global outputDir
    global threadLimit
    global pairendFlag
    global extsize
    global nofractionFlag
    global uniqueFlag
    global filterCutoff
    global shrink
    global prefix
    global padjMethod
    global testMethod
    ## report error log or not
    if finalFlag is True:
        errorFlag = True
    else:
        errorFlag = False
    ## get non-peak region row list
    nonPeakRegionRowList = FindNonPeakRegion(poolPeakRowList, geneAnnoFile, fragmentSize, tempDir)
    ## convert region row list into the SAF-format
    safMatrix = RegionToSafMatrix(poolPeakRowList, nonPeakRegionRowList, prefix + '.saf.matrix')
    ## construct final peak-reads dictionary
    bamFiles = ' '.join(bamFileList)
    featureCountFile = prefix + '.featureCounts.matrix'
    peakReadDict = RunFeatureCounts(safMatrix, featureCountFile, bamFiles, threadLimit, libraryStrand, pairendFlag, extsize, nofractionFlag, uniqueFlag, errorFlag)
    ## construct the sample matrix requried by DESeq2Peak.R
    sampleMtx = prefix + '.sample.matrix'
    with open(sampleMtx, 'w') as temp:
        row = ['sample', 'antibody\n']
        temp.write('\t'.join(row))
        for uniqBamIndex in bamDict['indexToInfo'].keys():
            ipType = bamDict['indexToInfo'][uniqBamIndex]['ipType']
            row = [uniqBamIndex, ipType + '\n']
            temp.write('\t'.join(row))
    ## construct the peak-counts matrix
    peakCountsMtx = prefix + '.counts.matrix'
    with open(peakCountsMtx, 'w') as temp:
        row = ['gene_id'] + uniqBamIndexList
        temp.write('\t'.join(row) + '\n')
        for peakid in sorted(peakReadDict.keys()):
            row = [peakid]
            for uniqBamIndex in uniqBamIndexList:
                counts = peakReadDict[peakid][uniqBamIndex]
                row.append(str(counts))
            temp.write('\t'.join(row) + '\n')
    ## run DESeq2Peak.R, prefix.DESeq2.txt, prefix.normalized.counts.txt
    if finalFlag is False:
        plotArg = '--noplot'
        runShrink = 'none'
    else:
        plotArg = ''
        runShrink = shrink
    deseq2Prefix = os.path.basename(prefix)
    if testMethod == 'Wald':
        command = 'DESeq2Peak.R --sampleMtx {} --counts {} \
            --control "input" --treat "ip" --design "antibody" --mean {} --norcounts \
            --shrink "{}" --prefix "{}" --relevel "antibody:input" \
            --padjust {} --test {} --output {} {}'
    else:
        command = 'DESeq2Peak.R --sampleMtx {} --counts {} \
            --control "input" --treat "ip" --design "antibody" --mean {} --norcounts \
            -shrink "{}" --prefix "{}" --relevel "antibody:input" --reduce "1" \
            --padjust {} --test {} --output {} {}'
    command = command.format(sampleMtx, peakCountsMtx, filterCutoff, runShrink, deseq2Prefix, padjMethod, testMethod, outputDir, plotArg)
    __ = SubprocessToList(command, errorFlag)
    ## result file used for downstream analysis
    norcountsFile = prefix + '.normalized.counts.txt'
    finalResFile = prefix + '.DESeq2.txt'
    ## delete temporary files
    if finalFlag is True and keepTempFlag is False:
        ## temporary files
        tempFileList = [safMatrix, sampleMtx, peakCountsMtx, norcountsFile, finalResFile]
        tempFileList += [prefix + '.MA.pdf', prefix + '.peak.MA.pdf', prefix + '.pvalue.histogram.pdf']
        tempFileList += [prefix + '.peak.pvalue.histogram.pdf', prefix + '.pvalueNorCounts.bar.pdf']
        # delete temporary file
        for tempFile in tempFileList:
            DeleteFile(tempFile)
    ## return results file
    return [finalResFile, norcountsFile]

########### main program ####################
## public viariables
outputDir = os.path.realpath(args.output)
fractionCutoff = args.fraction
minLength = args.min
maxLength = args.max
geneAnnoFile = args.anno
threadLimit = args.thread
pairendFlag = args.pairend
extsize = args.extsize
nofractionFlag = args.nofraction
uniqueFlag = args.unique
keepTempFlag = args.keepTemp
filterCutoff = args.filter
shrink = args.shrink
prefix = os.path.join(outputDir, args.prefix)
padjMethod = args.padjMethod
testMethod = args.test

if args.tempDir is not None:
    tempDir = os.path.realpath(args.tempDir)
else:
    tempDir = outputDir

## delete any tmp files generated by tempfile
tmpFileList = sorted(glob(os.path.join(outputDir, '*.tmp')) + glob(os.path.join(outputDir, 'sort*')))
if len(tmpFileList) > 0:
    for tmpFile in tmpFileList:
        DeleteFile(tmpFile)
## the number means the priority of software,
## in general, exomePeak and MeTPeak seems to be more accurate than others, so need to keep their tx structure first
poolCombination = 'poolCombination'
softwareDict = {'poolCombination':0, 'exomePeak':1, 'MeTPeak':1, 'exomePeak2':2, 'MACS2':3, 'SICER2':4, 'MoAIMS':5}

if len(args.name) != len(args.peak) or len(args.peak) != len(args.software):
    print('The length of --peak, --name and --software should be the same!')
    sys.exit()
else:
    for software in args.software:
        if software not in softwareDict:
            print('Only these programs ({}) are supported!'.format(','.join(sorted(softwareDict.keys()))))
            sys.exit()

bamFlagBool = True
if args.ip is not None and args.input is not None:
    if (len(set(args.ip)) > 1 and len(set(args.input)) == 0) or (len(set(args.input)) > 1 and len(set(args.ip)) == 0):
        print('The length of --ip or --input should be larger than 2!')
        sys.exit()
else:
    bamFlagBool = False

if bamFlagBool is True:
    if args.anno is None:
        print('--anno is required by --ip and --input!')
        sys.exit()
    else:
        if os.path.exists(args.anno) is False:
            print('--anno is required by --ip and --input!')
            sys.exit()
    args.ip = list(set(args.ip))
    args.input = list(set(args.input))

ipTypeList  = ['ip', 'input']
if bamFlagBool is True:
    ##build bam dict
    bamList = [args.ip, args.input]
    uniqBamIndexList = []
    bamDict = defaultdict(dict)
    bamDict['bamToIndex'] = defaultdict(dict)
    bamDict['indexToInfo'] = defaultdict(dict)
    for i in range(2):
        ipType = ipTypeList[i]
        bamFileList = bamList[i]
        for j in range(len(bamFileList)):
            bamFile = bamFileList[j]
            uniqBamIndex = '_'.join([ipType, str(j)])
            bamDict['bamToIndex'][bamFile] = uniqBamIndex
            bamDict['indexToInfo'][uniqBamIndex]['ipType'] = ipType
            uniqBamIndexList.append(uniqBamIndex)

# sort peak file by software priority
peakNameSep = '='
peakFileLabelDict = defaultdict(dict)
for i in range(len(args.peak)):
    bed = args.peak[i]
    name = args.name[i]
    software = args.software[i]
    peakFileLabelDict[bed] = [name, software]

peakFileList = sorted(peakFileLabelDict.keys(), key=lambda x:softwareDict[peakFileLabelDict[x][1]])
peakNameList = list(map(lambda x:peakFileLabelDict[x][0], peakFileList))
peakSoftwareList = list(map(lambda x:peakFileLabelDict[x][1], peakFileList))

# construct bedPeakRowList
## peakRow = [chr, start, end, ... , blockStarts] or [chr, start, end, name, score, strand]
## peakRowList =[bedrow, bedrow,...]
## bedPeakRowList = [peakRowList, peakRowList, ...]
bedPeakRowList = []
regexSkip = re.compile(r'^#')
for i in range(len(peakFileList)):
    bed = peakFileList[i]
    peakRowList = []
    with open(bed, 'r') as f:
        for line in f:
            if bool(regexSkip.match(line)):
                continue
            row = line.strip().split('\t')[0:12]
            peakBed = bedutils.buildbed(row)
            ## filter out peaks that larger or smaller than required
            if peakBed.exonlength > args.max or peakBed.exonlength < args.min:
                continue
            peakRowList.append(row)
    bedPeakRowList.append(peakRowList)

# construct bedLabelRowList
## labelRow = [rep1, exomePeak]
## labelRowList =[labelRow, labelRow,...]
## bedLabelRowList = [labelRowList, labelRowList, ...]
bedLabelRowList = []
for i in range(len(bedPeakRowList)):
    labelRowList = [[peakNameList[i], peakSoftwareList[i]] for k in bedPeakRowList[i]]
    bedLabelRowList.append(labelRowList)

## get intersected peaks between any two input peaks
indexList = list(range(len(bedPeakRowList)))
##[0,1,2,3] to [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
indexCombinationList = sorted(list(map(lambda x:sorted(list(x)), combinations(indexList, 2))), key=lambda x:[x[0], x[1]])

resultList = []
pool1 = Pool(processes=args.thread)
for index in range(len(indexCombinationList)):
    i,j = indexCombinationList[index]
    peak1RowList = bedPeakRowList[i]
    label1RowList = bedLabelRowList[i]
    peak2RowList = bedPeakRowList[j]
    label2RowList = bedLabelRowList[j]
    #peakRowList, labelRowList = PoolPeak(peakNameSep, peak1RowList, peak2RowList, label1RowList, label2RowList, softwareDict, args.min, args.max, tempDir, False)
    result = pool1.apply_async(PoolPeak, args=(peak1RowList, peak2RowList, label1RowList, label2RowList, False))
    resultList.append([index, result])
pool1.close()
pool1.join()

randomRestList = []
pool2 = Pool(processes=args.thread)
for result in resultList:
    index, poolResult = result
    peakRowList, labelRowList = poolResult.get()
    keyword = 'combinations-{}'.format(index)
    result = pool2.apply_async(RandomPool, args=(keyword, peakRowList, labelRowList, args.random))
    randomRestList.append([index, result])
pool2.close()
pool2.join()

comPeakRowList = [[] for i in indexCombinationList]
comLabelRowList = [[] for i in indexCombinationList]
for result in randomRestList:
    index, poolResult = result
    peakRowList, labelRowList = poolResult.get()
    comPeakRowList[index] = peakRowList
    comLabelRowList[index] = labelRowList

# get intersected comPeakRowList and comLabelRowList
comPeakRowList, comLabelRowList = CollapseDupBedRow(comPeakRowList, comLabelRowList, peakNameSep)
# continue to intersect peaks from peaks, to get pool peaks
peak1RowList = comPeakRowList[0]
label1RowList = comLabelRowList[0]
for i in range(1, len(comPeakRowList)):
    peak2RowList = comPeakRowList[i]
    label2RowList = comLabelRowList[i]
    peak1RowList, label1RowList = PoolPeak(peak1RowList, peak2RowList, label1RowList, label2RowList, False)

del comPeakRowList
del comLabelRowList

# collapse any possible overlapped peaks in the peak1RowList
peak1RowList, label1RowList = RandomPool('combination pool', peak1RowList, label1RowList, args.random)

# final round with original peaks
# append poolCombination to each label software, priority of poolCombination will set to 0
for i in range(len(label1RowList)):
    softwareList = label1RowList[i][1].split(',')
    softwareList.append(poolCombination)
    software = ','.join(softwareList)
    label1RowList[i][1] = software
# run the pool peaks with original peaks again to find any missed peaks by previous steps
for i in range(0, len(bedPeakRowList)):
    peak2RowList = bedPeakRowList[i]
    label2RowList = bedLabelRowList[i]
    peak1RowList, label1RowList = PoolPeak(peak1RowList, peak2RowList, label1RowList, label2RowList, True)
# remove poolCombination from each label software
for i in range(len(label1RowList)):
    software = ','.join(filter(lambda x:x != poolCombination, label1RowList[i][1].split(',')))
    label1RowList[i][1] = software
# collapse any possible overlapped peaks in the peak1RowList
keyword = 'pool with original peaks (final)'.format(i+1)
peak1RowList, label1RowList = RandomPool(keyword, peak1RowList, label1RowList, args.random)

peakLabelDict = {}
poolPeakRowList = []
peakLengthList = []
count = 1
for i in range(len(peak1RowList)):
    peakRow = peak1RowList[i]
    ## get exon length
    peakBed = bedutils.buildbed(peakRow)
    ## filter peaks larger than --max and --min
    if peakBed.exonlength > args.max or peakBed.exonlength < args.min:
        continue
    peakRow[3] = peakNameSep.join([peakRow[3], str(i + 1)])
    ## make thickStart == start and thickEnd == end
    peakRow[6] = peakRow[1]
    peakRow[7] = peakRow[2]
    ## get rid of poolCombination from labelRow
    label1RowList[i][1] = ','.join(filter(lambda x:x != poolCombination, label1RowList[i][1].split(',')))
    peakLabelDict[peakRow[3]] = label1RowList[i]
    poolPeakRowList.append(peakRow)
    peakLengthList.append(peakBed.exonlength)
## calculate average length of peak length
#peakAveLength = statistics.mean(peakLengthList)
# prepare peak dictionary
peakDict = defaultdict(dict)
for i in range(len(poolPeakRowList)):
    peakRow = poolPeakRowList[i]
    peakid = peakRow[3]
    peakDict[peakid]['row'] = peakRow
    peakBed = bedutils.buildbed(peakRow)
    peakDict[peakid]['length'] = peakBed.exonlength

## calculate reads from samples that cover peaks
if args.library == 'unstranded':
    libraryStrand = '0'
elif args.library == 'fr-firstrand':
    libraryStrand = '2'
else:
    libraryStrand = '1'

if bamFlagBool is True:
    bamFileList = args.ip + args.input
    ## get fragment length
    if args.pairend is False:
        readLengthList = []
        for bamFile in bamFileList:
            readLength = GetReadLength(bamFile)
            readLengthList.append(readLength)
        fragmentSize = max(readLengthList)
    else:
        fragmentSize = args.fragmentSize
    ## run DESeq2, first round: to keep significant peaks with log2fc larger than 0.05
    ## to remove the disturbance of non-peak region (background is not true)
    Deseq2ResFile, __ = RunDESeq2(poolPeakRowList, fragmentSize, bamFileList, libraryStrand, False)
    peakidSigFlagDict = {}
    with open(Deseq2ResFile, 'r') as f:
        row = f.readline().strip().split('\t')
        row = ['uniqId'] + row
        headerDict = defaultdict(dict)
        for i in range(len(row)):
            headerDict[row[i]] = i
        for line in f:
            row = line.strip().split('\t')
            peakid = row[headerDict['uniqId']]
            if peakid not in peakDict:
                continue
            log2fc = row[headerDict['log2FoldChange']]
            pval = row[headerDict['pvalue']]
            if pval != 'NA' and log2fc != 'NA':
                if float(log2fc) > 0:
                    peakidSigFlagDict[peakid] = True
                else:
                    peakidSigFlagDict[peakid] = False
            else:
                peakidSigFlagDict[peakid] = False
    ## get significant poolPeakRowList
    finalPoolPeakRowList = list(filter(lambda x: peakidSigFlagDict[x[3]] is True, poolPeakRowList))
    ## run DESeq2, to accuarately estimate the significant peaks with approriate background regions 
    Deseq2ResFile, norcountsFile = RunDESeq2(finalPoolPeakRowList, fragmentSize, bamFileList, libraryStrand, True)
    ## record peak-stats information
    peakStatsDict = defaultdict(dict)
    ## get normalized reads from prefix.normalized.counts.txt
    with open(norcountsFile, 'r') as f:
        ##ip_0, ip_1, input_0, input_1
        row = ['peakid'] + f.readline().strip().split('\t')
        rowLength = len(row)
        ipTypeIndexDict = {}
        for i in range(1, rowLength):
            uniqBamIndex = row[i]
            ipType = bamDict['indexToInfo'][uniqBamIndex]['ipType']
            ipTypeIndexDict[i] = ipType
        for line in f:
            row = line.strip().split('\t')
            peakid = row[0]
            if peakid not in peakDict:
                continue
            ipCountsDict = defaultdict(list)
            for i in range(1, rowLength):
                ipType = ipTypeIndexDict[i]
                ipCountsDict[ipType].append(float(row[i]))
            ## record average normalized counts
            peakDict[peakid]['ip'] = statistics.mean(ipCountsDict['ip'])
            peakDict[peakid]['input'] = statistics.mean(ipCountsDict['input'])
    ## record peak-log2fc information from prefix.DESeq2.txt file
    with open(Deseq2ResFile, 'r') as f:
        row = f.readline().strip().split('\t')
        row = ['uniqId'] + row
        headerDict = defaultdict(dict)
        for i in range(len(row)):
            headerDict[row[i]] = i
        for line in f:
            row = line.strip().split('\t')
            peakid = row[headerDict['uniqId']]
            if peakid not in peakDict:
                continue
            log2fc = row[headerDict['log2FoldChange']]
            pval = row[headerDict['pvalue']]
            pvalAdj = row[headerDict['padj']]
            ipMeanCounts = peakDict[peakid]['ip']
            inputMeanCounts = peakDict[peakid]['input']
            valueRow = [ipMeanCounts, inputMeanCounts, log2fc, pval, pvalAdj]
            peakStatsDict[peakid]['all'] = valueRow
            ## discard the peaks with log2fc <= 0 or pvalue larger than --paval
            peakStatsDict[peakid]['sigfc'] = False
            ## calculate ip density (number of ip reads per kb)
            #ipDensity = ipMeanCounts / peakDict[peakid]['length'] * 1000
            if pval != 'NA' and pvalAdj != 'NA' and log2fc != 'NA':
                if ipMeanCounts >= args.filter and float(pval) <= args.pval and float(pvalAdj) <= args.padj and float(log2fc) >= args.log2fc:
                    peakStatsDict[peakid]['sigfc'] = True
else:
    finalPoolPeakRowList = poolPeakRowList

## generate final output results
headerRow = ['#chr', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
headerRow += ['thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
headerRow += ['peakLength', 'source', 'software']

if bamFlagBool is True:
    headerRow += ['mean.normalized.ip.counts', 'mean.normalized.input.counts']
    headerRow += ['log2FoldChange', 'pval.' + args.test, 'padj.' + args.padjMethod, 'ipDensity']

outputRowList = []
## peakid -> row
for i in range(len(finalPoolPeakRowList)):
    peakRow = finalPoolPeakRowList[i]
    peakid = peakRow[3]
    peakRow[3] = peakid.split(peakNameSep)[0]
    peakLength = peakDict[peakid]['length']
    labelList = list(map(lambda x: ','.join(sorted(set(x.split(',')))), peakLabelDict[peakid]))
    peakRow.append(peakLength)
    peakRow.extend(labelList)
    if bamFlagBool is True:
        if peakid not in peakStatsDict:
            continue
        row[4] = peakStatsDict[peakid]['all'][3]
        if row[4] == 'NA':
            row[4] = 1
        row += peakStatsDict[peakid]['all']
    outputRowList.append([peakid, row])

## output results
allOutput = prefix + '.bed'
with open(allOutput, 'w') as out:
    out.write('\t'.join(map(str, headerRow)) + '\n')
    for outputRow in outputRowList:
        peakid = outputRow[0]
        row = outputRow[1]
        out.write('\t'.join(map(str, row)) + '\n')

## output significantly changed results
sigfcOutput = prefix + '.sigfc.bed'
if bamFlagBool is True:
    with open(sigfcOutput, 'w') as out:
        out.write('\t'.join(map(str, headerRow)) + '\n')
        for outputRow in outputRowList:
            peakid = outputRow[0]
            if peakStatsDict[peakid]['sigfc'] is True:
                row = outputRow[1]
                out.write('\t'.join(map(str, row)) + '\n')
