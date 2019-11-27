#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
import tempfile
import subprocess
## own module
import BedMan

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-anno', action='store', type=str, required=True,
                    help='Gene annotation file in bed12 format (main annotation) \
                    (4th column:gene_id:gene_name:gene_type:tx_id:tx_name:tx_type)')
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input peak file (bed6 or bed6+)')
parser.add_argument('-bf', action='store', type=str,
                    help='bedtools: -f')
parser.add_argument('-bF', action='store', type=str,
                    help='bedtools: -F')
parser.add_argument('-br', action='store', type=str,
                    help='bedtools: -r')
parser.add_argument('-be', action='store', type=str,
                    help='bedtools: -e')
parser.add_argument('-method', action='store', type=str, choices=['center', 'border'],
                    default='border',
                    help='use peak center|border to calcualte distance to TSS|TTS')
parser.add_argument('-name', action='store', type=str,
                    default='peak=',
                    help='prefix for each peak name')
parser.add_argument('-kbTSS', action='store', type=str,
                    default='1,2,3,5:1,2',
                    help='<int,int,...:int,int,...> record (up:down) distance to TSS in kb')
parser.add_argument('-kbTTS', action='store', type=str,
                    default='1:1,2',
                    help='<int,int,...:int,int,...> record (up:down) distance to TSS in kb')
parser.add_argument('-codon', action='store', type=str,
                    help='<int#1,int#2>#1 bp upstream and #2 bp downstream start/stop codon as their feature')
parser.add_argument('-mode', action='store', type=str, choices=['DNA', 'RNA'],
                    default='RNA',
                    help='Annotation mode')
parser.add_argument('-strand', action='store_true',
                    default=False,
                    help='bedtools: -s (strand)')
parser.add_argument('--disablePeakUniq', action='store_true',
                    default=False,
                    help='output multiple records for same peak')
parser.add_argument('--keepName', action='store_true',
                    default=False,
                    help='keep original 4th column in output')
parser.add_argument('-extraAnno', nargs='*', type=str,
                    help='individual extraitional annotations in bed6 format \
                    (4th column:identifier1:identifier2:identifier3:type) \
                    (the order represents type priority)')
parser.add_argument('-extraType', nargs='*', type=str, choices=['gene', 'element'],
                    help='type of extraitional annotations (gene|element, same length with -extraAnno)')
parser.add_argument('-gsize', action='store', type=str,
                    help='genome size file (required for DNA mode)')
parser.add_argument('-priF', action='store', type=str,
                    help="Defined priority of features (eg., \"CDS,5' UTR\")")
parser.add_argument('-priG', action='store', type=str,
                    help='Defined priority of gene types (eg., "protein-coding,lncRNA")')
parser.add_argument('-geneClassFile', action='store', type=str, required=True,
                    help='geneType-geneClass pairwise file')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output result (coordiantes in 0-base)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

# public arguments
mainAnno = args.anno
peak = args.input
output = args.output
uniqPeakPre = args.name

# decode codon, kbTSS, kbTTS region
codonExtList = list()
kbTssList = list()
kbTtsList = list()
if bool(args.codon):
    try:
        codonExtList = list(map(int, args.codon.split(',')))
        if len(codonExtList) != 2:
            sys.stderr.write('Incorrect -codon!')
            sys.exit()
    except ValueError as e:
        sys.stderr.write('Incorrect -codon!')
        sys.exit()

try:
    kbTssList = list(map(lambda x: [int(y) for y in x], 
        [z.split(',') for z in args.kbTSS.split(':')]))
    if len(kbTssList) != 2:
        sys.stderr.write('Incorrect -kbTSS!')
        sys.exit()
except ValueError as e:
    sys.stderr.write('Incorrect -kbTSS!')
    sys.exit()

try:
    kbTtsList = list(map(lambda x: [int(y) for y in x], 
        [z.split(',') for z in args.kbTTS.split(':')]))
    if len(kbTtsList) != 2:
        sys.stderr.write('Incorrect -kbTTS!')
        sys.exit()
except ValueError as e:
    sys.stderr.write('Incorrect -kbTTS!')
    sys.exit()

# defined priority of features
if bool(args.priF):
    priFeatureList = args.priF
else:
    if args.mode == 'RNA':
        priFeatureList = ['CDS', 'start_codon', 'stop_codon', "5' UTR", "3' UTR", 'Exon', 'Intron']
    else:
        priFeatureList = ["5' UTR", 'start_codon', 'CDS', 'stop_codon', "3' UTR", 'Exon', 'Intron']
        mainKbTssList = [0] + sorted(set(kbTssList[0] + kbTssList[1]))
        priTssList = list()
        # format TSS (<=1kb), TSS (1-2kb)
        for i in range(len(mainKbTssList)-1):
            if i == 0:
                feature = 'TSS (<={0}kb)'.format(mainKbTssList[1])
            else:
                feature = 'TSS ({0}-{1}kb)'.format(mainKbTssList[i], mainKbTssList[i+1])
            priTssList.append(feature)
        mainKbTtsList = [0] + sorted(set(kbTtsList[0] + kbTtsList[1]))
        priTtsList = list()
        for i in range(len(mainKbTtsList)-1):
            if i == 0:
                feature = 'TTS (<={0}kb)'.format(mainKbTtsList[1])
            else:
                feature = 'TTS ({0}-{1}kb)'.format(mainKbTtsList[i], mainKbTtsList[i+1])
            priTtsList.append(feature)
        priFeatureList = priTssList + priFeatureList + priTtsList
        priFeatureList.append('intergenic')

priFeatureNum = len(priFeatureList)

if bool(args.priG):
    priGenetList = args.priG
else:
    priGenetList = ['protein_coding', "lncRNA", "TR_gene", 'IG_gene', 'miRNA',\
        'snoRNA', 'rRNA', 'sncRNA', 'pseudogene', 'other']
priGenetNum = len(priGenetList)

# annotate peak in DNA|RNA mode
def annoPeak(mode, peakBed, annoBed, bf, bF, br, be, bs, annoType):
    if mode == 'RNA':
        command = 'bedtools intersect -a {peakBed} -b {annoBed} \
            -wa -wb {bf} {bF} {br} {be} {bs}'.format(**vars())
        annoResList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
    else:
        if annoType == 'element':
            command = 'bedtools intersect -a {peakBed} -b {annoBed} \
                -wa -wb {bf} {bF} {br} {be}'.format(**vars())
            annoResList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
        else:
            annoExtBedTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
            ## extract TSS
            awkCommand = 'awk \'BEGIN{OFS="\t";FS="\t"}{if($6=="+"){\
                print($1,$2,$2+1,$4":TSS:"$2,$5,$6)}else{print($1,$3-1,$3,$4":TSS:"$3-1,$5,$6)}}\' ' + annoBed
            ## slop TSS
            farmostUpTss = kbTssList[0][-1]*1000
            farmostDownTss = kbTssList[1][-1]*1000
            slopCommand = 'bedtools slop -i stdin -l {farmostUpTss} -r {farmostDownTss} -s -g '.format(**vars()) + args.gsize
            intersectComand = 'bedtools intersect -a {peakBed} -b stdin -wa -wb {bf} {bF} {br} {be}'.format(**vars())
            dnaTssAnnoCommand = "{awkCommand} | {slopCommand} | {intersectComand}".format(**vars())
            dnaTssAnnoResList = bytes.decode(subprocess.check_output(dnaTssAnnoCommand, shell=True)).split('\n')
            ## extract TTS
            awkCommand = 'awk \'BEGIN{OFS="\t";FS="\t"}{if($6=="+"){\
                print($1,$3-1,$3,$4":TTS:"$3-1,$5,$6)}else{print($1,$2,$2+1,$4":TTS:"$2,$5,$6)}}\' ' + annoBed
            ## slop TTS
            farmostUpTts = kbTtsList[0][-1]
            farmostDownTts = kbTtsList[1][-1]
            slopCommand = 'bedtools slop -i stdin -l {farmostUpTts} -r {farmostDownTts} -s -g '.format(**vars()) + args.gsize
            intersectComand = 'bedtools intersect -a {peakBed} -b stdin -wa -wb {bf} {bF} {br} {be}'.format(**vars())
            dnaTtsAnnoCommand = "{awkCommand} | {slopCommand} | {intersectComand}".format(**vars())
            dnaTtsAnnoResList = bytes.decode(subprocess.check_output(dnaTtsAnnoCommand, shell=True)).split('\n')
            annoResList = dnaTssAnnoResList + dnaTtsAnnoResList
    return annoResList

# decode functions
def rnaFeatureDecode (bed12Row, peakLocus, distBool=False):
    ## decode bed12 with BedMan.decodeBed12, dispite strand
    ## other-RNA:[[exonL], [intronL]]
    ## mRNA:[[exonL], [intronL], [[thickUpL], [thickInL], [thickDownL]]]
    strand = bed12Row[5]
    decodeList = BedMan.decodeBed12(bed12Row)
    annoFeatureDict = defaultdict(list)
    priorityIndexList = list()
    if len(decodeList) == 2:
        exonL, intronL = decodeList
        ### reverse to get true number if strand '-'
        if strand == '-':
            exonL.reverse()
            intronL.reverse()
        featureDict = defaultdict(list)
        ### label overlap exon|intron, using BedMan.overlap
        featureDict['Exon'] = list(map(lambda x:BedMan.overlap(peakLocus, x), exonL))
        featureDict['Intron'] = list(map(lambda x:BedMan.overlap(peakLocus, x), intronL))
    else:
        exonL, intronL, thickL = decodeList
        ### reverse to get true number if strand '-'
        if strand == '-':
            intronL.reverse()
            ### reverse thickL to get 5' UTR, cds, 3' UTR
            thickL.reverse()
            for i in range(len(thickL)):
                thickL[i].reverse()
        utr5L, cdsL, utr3L = thickL
        ### label overlap 5'utr|cds|3'utr|intron, using BedMan.overlap
        featureDict = defaultdict(list)
        if bool(utr5L):
            featureDict["5' UTR"] = list(map(lambda x:BedMan.overlap(peakLocus, x), utr5L))
        if bool(utr3L):
            featureDict["3' UTR"] = list(map(lambda x:BedMan.overlap(peakLocus, x), utr3L))
        featureDict['CDS'] = list(map(lambda x:BedMan.overlap(peakLocus, x), cdsL))
        featureDict['Intron'] = list(map(lambda x:BedMan.overlap(peakLocus, x), intronL))
        ### if args.codon, defined start&stop codon feature
        if bool(codonExtList):
            if strand == '+':
                startCodonL = [[int(bed12Row[6]) - codonExtList[0], int(bed12Row[6]) + codonExtList[1]]]
                stopCodonL = [[int(bed12Row[7]) - codonExtList[0], int(bed12Row[7]) + codonExtList[1]]]
            else:
                stopCodonL = [[int(bed12Row[6]) - codonExtList[0], int(bed12Row[6]) + codonExtList[1]]]
                startCodonL = [[int(bed12Row[7]) - codonExtList[0], int(bed12Row[7]) + codonExtList[1]]]
            featureDict['start_codon'] = list(map(lambda x:BedMan.overlap(peakLocus, x), startCodonL))
            featureDict['stop_codon'] = list(map(lambda x:BedMan.overlap(peakLocus, x), stopCodonL))
    ## index priority
    ## sort by order in priFeatureList 
    featureList = sorted(featureDict.keys(), key = lambda x:priFeatureList.index(x))
    ## filter out feature with no overlaps
    featureList = list(filter(lambda x: sum(featureDict[x]), featureList))
    priorityIndexList = list(map(lambda x:priFeatureList.index(x), featureList))
    ## mRNA:[1, 0, 2, 4, 3]
    minIndex = priorityIndexList.index(min(priorityIndexList))
    annoFeatureDict['main'] = featureList[minIndex]
    annoFeatureDict['minor'] = '+'.join(featureList)
    detailList = list()
    for minor in featureList:
        orderList = list()
        for i in range(len(featureDict[minor])):
            if featureDict[minor][i] > 0:
                orderList.append(i+1)
        annoFeatureDict['detail'].extend(list(map(lambda x:minor+'-'+str(x), orderList)))
    ## calculate distance to TSS for DNA mode
    if distBool:
        tss = bed12Row[1]
        if strand == '-':
            tss = bed12Row[2]
        tss = int(tss)
        peakCoor = peakLocus[0]
        if args.method == 'center':
            peakCoor = int(sum(peakLocus) / 2)
        distance = peakCoor - tss
        if strand == '-':
            distance = - distance
        annoFeatureDict['distance'] = distance
    ## return annoated region information
    return annoFeatureDict

# decode DNA fetures
def dnaFeatureDecode (bedAnnoRow, peakLocus, annoType):
    annoFeatureDict = defaultdict(dict)
    strand = bedAnnoRow[5]
    txInfo = bedAnnoRow[3].split(':')
    if annoType == 'element':
        annoFeatureDict['main'] = txInfo[-1]
        annoStart = int(bedAnnoRow[1])
        annoEnd = int(bedAnnoRow[2])
        if strand == '.':
            ## calcualte distance, return min(overlap, length-overlap)
            annoLength = annoEnd - annoStart
            overlapLength = BedMan.overlap([annoStart, annoEnd], peakLocus)
            distance = min(overlapLength, (annoLength - overlapLength))
        else:
            featureCoor = annoStart
            if strand == '-':
                featureCoor = annoEnd
            ## calulating distance
            if args.method == 'border':
                if featureCoor >= peakLocus[1]:
                    distance = peakLocus[1] - 1 - featureCoor
                elif featureCoor <= peakLocus[0]:
                    distance = peakLocus[0] - featureCoor
                else:
                    if abs(peakLocus[0] - featureCoor) >= (peakLocus[1] - 1 - featureCoor):
                        distance = peakLocus[0] - featureCoor
                    else:
                        distance = peakLocus[1] - 1 - featureCoor
            else:
                centerPeakCoor = int(sum(peakLocus) / 2)
                distance = centerPeakCoor - featureCoor
            if strand == '-':
                distance = - distance
        annoFeatureDict['distance'] = distance
    else:
        ## TSS or TTS
        featureCoor = int(txInfo[-1])
        mainFeature = txInfo[-2]
        ## calulating distance
        if args.method == 'border':
            if featureCoor >= peakLocus[1]:
                distance = peakLocus[1] - 1 - featureCoor
            elif featureCoor <= peakLocus[0]:
                distance = peakLocus[0] - featureCoor
            else:
                if abs(peakLocus[0] - featureCoor) >= (peakLocus[1] - 1 - featureCoor):
                    distance = peakLocus[0] - featureCoor
                else:
                    distance = peakLocus[1] - 1 - featureCoor
        else:
            centerPeakCoor = int(sum(peakLocus) / 2)
            distance = centerPeakCoor - featureCoor
        if strand == '-':
            distance = - distance
        ## determin feature ground on distance
        if mainFeature == 'TSS':
            mainKbList = [0] + sorted(set(kbTssList[0] + kbTssList[1]))
        else:
            mainKbList = [0] + sorted(set(kbTtsList[0] + kbTtsList[1]))
        if abs(distance) <= mainKbList[-1]*1000:
            for i in range(len(mainKbList)-1):
                if abs(distance) >= mainKbList[i]*1000 and abs(distance) <= mainKbList[i+1]*1000:
                    if i == 0:
                        feature = '{0} (<={1}kb)'.format(mainFeature, mainKbList[1])
                    else:
                        feature = '{0} ({1}-{2}kb)'.format(mainFeature, mainKbList[i], mainKbList[i+1])
                    break
            annoFeatureDict['main'] = feature
        else:
            annoFeatureDict['main'] = 'intergenic'
        annoFeatureDict['distance'] = distance
    return annoFeatureDict

# running main annotation
bs = '-s' if args.strand else ''
bf = '-f ' + args.bf if bool(args.bf) else ''
bF = '-F ' + args.bF if bool(args.bF) else ''
br = '-r ' + args.br if bool(args.br) else ''
be = '-e ' + args.be if bool(args.be) else ''

# construct bed6 from input
peakIdList = list()
count = 1
skipList = list()
peakInfoDict = defaultdict(dict)
peakBed6Tmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
peakHeaderRow = list()
peakColNum = 0
with open(args.input, 'r') as f, open(peakBed6Tmp.name, 'w') as temp:
    for line in f:
        ## extract peak header, if have
        if count == 1:
            if bool(re.match(r'^#', line)):
                peakHeaderRow = line.strip().split('\t')
                continue
        ## skip and report non-valid bed6 line
        if BedMan.check(line, ctype='bed6') is False:
            skipList.append(count)
            continue
        row = line.strip().split('\t')
        if peakColNum == 0:
            peakColNum = len(row)
        originalId = row[3]
        uniqPeakId = uniqPeakPre + str(count)
        peakIdList.append(uniqPeakId)
        ## keep original peak record
        peakInfoDict[uniqPeakId]['original'] = originalId
        peakInfoDict[uniqPeakId]['row'] = row
        ## construct peak-bed6
        row[3] = uniqPeakId
        bed6Line = '\t'.join(row[:6]) + '\n'
        temp.write(bed6Line)
        count += 1
peakBed6Tmp.seek(0)

# running main annotation
mainAnnoResList = annoPeak('RNA', peakBed6Tmp.name, mainAnno, bf, bF, br, be, bs, 'gene')

# DNA mode: generate promoter and enhancer feature
if args.mode == 'DNA':
    if bool(args.gsize) is False:
        sys.stderr.write('No -gsize!')
        sys.exit()
    dnaAnnoResList = annoPeak('DNA', peakBed6Tmp.name, mainAnno, bf, bF, br, be, bs, 'gene')
    mainAnnoResList = mainAnnoResList + dnaAnnoResList

# built up gene-type pairwise relationships
geneClassDict = defaultdict(dict)
with open(args.geneClassFile) as f:
    for line in f:
        row = line.strip().split('\t')
        mainType = row[0]
        geneType = row[1]
        geneClassDict[geneType] = mainType
# built up tx-infor pairwise relationships
txDict = defaultdict(dict)
with open(args.anno) as f:
    for line in f:
        row = line.strip().split('\t')
        txInfo = row[3].split(':')
        txId = txInfo[3]
        txDict[txId]['geneId'] = txInfo[0]
        txDict[txId]['geneName'] = txInfo[1]
        txDict[txId]['geneType'] = txInfo[2]
        txDict[txId]['txName'] = txInfo[4]
        txDict[txId]['txType'] = txInfo[5]
        txDict[txId]['strand'] = row[5]
        txDict[txId]['txStart'] = int(row[1])
        txDict[txId]['txEnd'] = int(row[2])

# decode peak-anno pairwise relationships
mainAnnoPeakDict = defaultdict(dict)
for line in mainAnnoResList:
    if bool(line) is False:
        continue
    row = line.strip().split('\t')
    peakRow = row[0:6]
    bedAnnoRow = row[6:]
    ## get peak locus
    peakRow[1] = int(peakRow[1])
    peakRow[2] = int(peakRow[2])
    peakId = peakRow[3]
    peakLocus = [peakRow[1], peakRow[2]]

    txInfo = bedAnnoRow[3].split(':')
    geneId = txInfo[0]
    geneName = txInfo[1]
    geneType = txInfo[2]
    txId = txInfo[3]
    if args.mode == 'RNA':
        annoFeatureDict = rnaFeatureDecode(bedAnnoRow, peakLocus, distBool=False)
    else:
        if len(row) == 12:
            annoFeatureDict = dnaFeatureDecode(bedAnnoRow, peakLocus, 'gene')
        else:
            annoFeatureDict = rnaFeatureDecode(bedAnnoRow, peakLocus, distBool=True)
    if peakId in mainAnnoPeakDict:
        if geneId in mainAnnoPeakDict[peakId]:
            mainAnnoPeakDict[peakId][geneId]['feature'].append(annoFeatureDict)
            mainAnnoPeakDict[peakId][geneId]['tx'].append(txId)
        else:
            mainAnnoPeakDict[peakId][geneId]['feature'] = [annoFeatureDict]
            mainAnnoPeakDict[peakId][geneId]['tx'] = [txId]
            mainAnnoPeakDict[peakId][geneId]['name'] = geneName
            mainAnnoPeakDict[peakId][geneId]['type'] = geneType
    else:
        mainAnnoPeakDict[peakId] = defaultdict(dict)
        mainAnnoPeakDict[peakId][geneId]['feature'] = [annoFeatureDict]
        mainAnnoPeakDict[peakId][geneId]['tx'] = [txId]
        mainAnnoPeakDict[peakId][geneId]['name'] = geneName
        mainAnnoPeakDict[peakId][geneId]['type'] = geneType

annoPeakDict = defaultdict(dict)
for peakId in mainAnnoPeakDict.keys():
    geneFeatureDict = defaultdict(dict)
    for geneId in mainAnnoPeakDict[peakId].keys():
        ## determin feature for peak-gene (multiple tx)
        annoFeatureDictList = mainAnnoPeakDict[peakId][geneId]['feature']
        txList = mainAnnoPeakDict[peakId][geneId]['tx']
        mainFeatureList = list(map(lambda x:x['main'], annoFeatureDictList))
        ## get index in priFeatureList
        featureIndexList = list(map(lambda x:priFeatureList.index(x), mainFeatureList))
        featureLabelList = list()
        ## get first priority index in priFeatureList
        minIndex = min(featureIndexList)
        for i in range(len(featureIndexList)):
            if featureIndexList[i] == minIndex:
                featureLabelList.append(i)
        label = featureLabelList[0]
        if len(featureLabelList) > 1:
            ## if multiple first priority, return the maximun length of tx
            if args.mode == 'RNA':
                txLengthList = list(map(lambda x: txDict[x]['txEnd'] - txDict[x]['txStart'], txList))
                txMaxLen = 0
                for i in featureLabelList:
                    if txLengthList[i] > txMaxLen:
                        txMaxLen = txLengthList[i]
                        label = i
            else:
                ## if multiple first priority, return the closest tx
                distList = list(map(lambda x: x['distance'], annoFeatureDictList))
                minDist = 10**8
                for i in featureLabelList:
                    if abs(distList[i]) < abs(minDist):
                        minDist = distList[i]
                        label = i
        ## get final tx information
        txId = txList[label]
        geneFeatureDict[geneId]['txId'] = txId
        txType = txDict[txId]['txType']
        geneFeatureDict[geneId]['txName'] = txDict[txId]['txName']
        geneFeatureDict[geneId]['txType'] = txType
        geneFeatureDict[geneId]['geneName'] = txDict[txId]['geneName']
        geneFeatureDict[geneId]['geneType'] = txDict[txId]['geneType']
        if txType in geneClassDict:
            geneFeatureDict[geneId]['geneClass'] = geneClassDict[txType]
        else:
            geneFeatureDict[geneId]['geneClass'] = 'other'
        ## get final feature
        geneFeatureDict[geneId]['main'] = annoFeatureDictList[label]['main']
        if args.mode == 'RNA':
            geneFeatureDict[geneId]['minor'] = annoFeatureDictList[label]['minor']
            geneFeatureDict[geneId]['detail'] = annoFeatureDictList[label]['detail']
        else:
            geneFeatureDict[geneId]['distance'] = annoFeatureDictList[label]['distance']
    ## store final annotation as dict
    if args.disablePeakUniq:
        annoPeakDict[peakId]['main'] = geneFeatureDict
    else:
        ## make meta-gene unique
        geneIdList = sorted(geneFeatureDict.keys())
        geneClassList = list(map(lambda x:geneFeatureDict[x]['geneClass'], geneIdList))
        geneClassIndexList = list(map(lambda x:priGenetList.index(x), geneClassList))
        if args.mode == 'RNA':
            ## if multiple first priority, return the maximun length of tx
            geneLabelList = list()
            ## get first priority index in priGenetList
            minIndex = min(geneClassIndexList)
            for i in range(len(geneClassIndexList)):
                if geneClassIndexList[i] == minIndex:
                    geneLabelList.append(i)
            label = geneLabelList[0]
            if len(geneLabelList) > 1:
                txList = list(map(lambda x:geneFeatureDict[x]['txId'], geneIdList))
                txLengthList = list(map(lambda x: txDict[x]['txEnd'] - txDict[x]['txStart'], txList))
                txMaxLen = 0
                for i in geneLabelList:
                    if txLengthList[i] > txMaxLen:
                        txMaxLen = txLengthList[i]
                        label = i
        else:
            ## if multiple first priority, return the closest tx
            distList = list(map(lambda x: geneFeatureDict[x]['distance'], geneIdList))
            minIndex = min(list(map(abs, distList)))
            distLabelList = list()
            for i in range(len(distList)):
                if abs(distList[i]) == minIndex:
                    distLabelList.append(i)
            label = distLabelList[0]
            priGenetypeIndex = 100
            if len(distLabelList) > 1:
                for i in distLabelList:
                    if geneClassIndexList[i] < priGenetypeIndex:
                        priGenetypeIndex = geneClassIndexList[i]
                        label = i
        ## get final main annotation
        labelGeneId = geneIdList[label]
        labelDict = defaultdict(dict)
        labelDict[labelGeneId] = geneFeatureDict[labelGeneId]
        annoPeakDict[peakId]['main'] = labelDict

## run extra annotation, if -extraition
if bool(args.extraAnno):
    extraAnnoList = args.extraAnno
    extraTypeList = args.extraType
    ## determin type by oder
    ## store annoPeak results
    extraPriTypeList = list()
    annoResList = list()
    for i in range(len(extraAnnoList)):
        extraAnno = extraAnnoList[i]
        extraType = extraTypeList[i]
        with open(extraAnno, 'r') as f:
            for line in f:
                if bool(re.match(r'^#', line)):
                    continue
                else:
                    row = line.strip().split('\t')
                    inforList = row[3].split(':')
                    extraAnnoType = inforList[3]
                    extraPriTypeList.append(extraAnnoType)
                    break
        eachAnnoResList = annoPeak(args.mode, peakBed6Tmp.name, extraAnno, bf, bF, br, be, bs, extraType)
        annoResList.append(eachAnnoResList)
    #extraAnnoPeakDict
    extraAnnoPeakDict = defaultdict(dict)
    for i in range(len(extraAnnoList)):
        eachAnnoResList = annoResList[i]
        extraType = extraTypeList[i]
        for line in eachAnnoResList:
            if bool(line) is False:
                continue
            try:
                row = line.strip().split('\t')
            except AttributeError as e:
                print(line)
            row = line.strip().split('\t')
            peakRow = row[0:6]
            bedAnnoRow = row[6:]
            ## get peak locus
            peakRow[1] = int(peakRow[1])
            peakRow[2] = int(peakRow[2])
            peakId = peakRow[3]
            peakLocus = [peakRow[1], peakRow[2]]
        
            infoList = bedAnnoRow[3].split(':')
            identifier1 = infoList[0]
            identifier2 = infoList[1]
            identifier3 = infoList[2]
            identifier = ':'.join([identifier1, identifier2, identifier3])
            extraAnnoType = infoList[3]
            ## decode feature
            if args.mode == 'RNA':
                annoFeatureDict = defaultdict(dict)
                annoFeatureDict['main'] = 'Exon'
            else:
                annoFeatureDict = dnaFeatureDecode(bedAnnoRow, peakLocus, extraType)
            extraAnnoPeakDict[peakId][identifier] = defaultdict(dict)
            extraAnnoPeakDict[peakId][identifier]['feature'] = annoFeatureDict
            extraAnnoPeakDict[peakId][identifier]['length'] = int(bedAnnoRow[2]) - int(bedAnnoRow[1])
            extraAnnoPeakDict[peakId][identifier]['extraType'] = extraAnnoType
    ## collapse feature and gene type
    for peakId in extraAnnoPeakDict.keys():
        if args.disablePeakUniq:
            annoPeakDict[peakId]['extra'] = extraAnnoPeakDict[peakId]
        else:
            tempDict = extraAnnoPeakDict[peakId]
            identifierList = sorted(tempDict.keys())
            typeList = list(map(lambda x:tempDict[x]['extraType'], identifierList))
            typeIndexList = list(map(lambda x:extraPriTypeList.index(x), typeList))
            if args.mode == 'RNA':
                ## get first priority index in extraPriTypeList
                typeLabelList = list()
                minIndex = min(typeIndexList)
                for i in range(len(typeIndexList)):
                    if typeIndexList[i] == minIndex:
                        typeLabelList.append(i)
                label = typeLabelList[0]
                if len(typeLabelList) > 1:
                    lengthList = list(map(lambda x: tempDict[x]['length'], identifierList))
                    maxLen = 0
                    for i in typeLabelList:
                        if lengthList[i] > maxLen:
                            maxLen = lengthList[i]
                            label = i
            else:
                ## if multiple first priority, return the closest tx
                distList = list(map(lambda x: tempDict[x]['feature']['distance'], identifierList))
                minIndex = min(list(map(abs, distList)))
                distLabelList = list()
                for i in range(len(distList)):
                    if abs(distList[i]) == minIndex:
                        distLabelList.append(i)
                label = distLabelList[0]
                priExtraAnnoTypeIndex = 100
                if len(distLabelList) > 1:
                    for i in distLabelList:
                        if typeIndexList[i] < priExtraAnnoTypeIndex:
                            priExtraAnnoTypeIndex = typeIndexList[i]
                            label = i
            labelId = identifierList[label]
            labelDict = defaultdict(dict)
            labelDict[labelId] = tempDict[labelId]
            annoPeakDict[peakId]['extra'] = labelDict

# ready for output
## construct header
headerRow = list()
if len(peakHeaderRow) != peakColNum:
    peakHeaderRow = ['peakCol'+str(i+1) for i in range(peakColNum)]
headerRow.extend(peakHeaderRow)
mainHeaderRow = ["GeneId", "GeneName", "GeneType", "GeneClass", "TxId", "TxName", "TxType", "Feature"]
extraHeaderRow = ['ExtraName-1', 'ExtraName-2', 'ExtraName-3', 'ExtraFeature', 'ExtraType']
if args.mode == 'RNA':
    mainHeaderRow.append('MinorFeature')
    mainHeaderRow.append('DetailFeature')
else:
    mainHeaderRow.append('Distance')
    extraHeaderRow.append('ExtraDistance')
headerRow.extend(mainHeaderRow)
headerRow.extend(extraHeaderRow)
## construct output contents
outputRowList = [headerRow]
for peakId in peakIdList:
    peakRow = peakInfoDict[peakId]['row']
    if args.keepName is False:
        peakRow[3] = peakId
    if peakId in annoPeakDict:
        extraAnnoRow = list()
        if 'extra' in annoPeakDict[peakId]:
            extraAnnoPeakDict = annoPeakDict[peakId]['extra']
            identifierList = sorted(extraAnnoPeakDict.keys())
            identifierRow = [[], [], []]
            featureRow = list()
            extraAnnoTypeRow = list()
            distanceRow = list()
            for identifier in sorted(extraAnnoPeakDict.keys()):
                identifier1, identifier2, identifier3 = identifier.split(':')[0:3]
                identifierRow[0].append(identifier1)
                identifierRow[1].append(identifier2)
                identifierRow[2].append(identifier3)
                featureRow.append(extraAnnoPeakDict[identifier]['feature']['main'])
                extraAnnoTypeRow.append(extraAnnoPeakDict[identifier]['extraType'])
                if args.mode == 'DNA':
                    distanceRow.append(extraAnnoPeakDict[identifier]['feature']['distance'])
            for i in range(len(identifierRow)):
                identifierRow[i] = ','.join(identifierRow[i])
            identifier1Col = identifierRow[0]
            identifier2Col = identifierRow[1]
            identifier3Col = identifierRow[2]
            try:
                featureCol = ','.join(extraAnnoPeakDict)
            except TypeError as e:
                print(featureRow)
            featureCol = ','.join(featureRow)
            extraAnnoTypeCol = ','.join(extraAnnoTypeRow)
            ## 5 or 6 elements
            extraAnnoRow = [identifier1Col, identifier2Col, identifier3Col, featureCol, extraAnnoTypeCol]
            if args.mode == 'DNA':
                distanceCol = ','.join(map(str,distanceRow))
                extraAnnoRow.append(distanceCol)
        else:
            if args.mode == 'RNA':
                extraAnnoRow = ['na' for i in range(5)]
            else:
                extraAnnoRow = ['na' for i in range(6)]
        mainAnnoRow = list()
        if 'main' in annoPeakDict[peakId]:
            mainAnnoPeakDict = annoPeakDict[peakId]['main']
            for geneId in sorted(mainAnnoPeakDict.keys()):
                geneName = mainAnnoPeakDict[geneId]['geneName']
                geneType = mainAnnoPeakDict[geneId]['geneType']
                geneClass = mainAnnoPeakDict[geneId]['geneClass']
                txId = mainAnnoPeakDict[geneId]['txId']
                txName = mainAnnoPeakDict[geneId]['txName']
                txType = mainAnnoPeakDict[geneId]['txType']
                mainFeature = mainAnnoPeakDict[geneId]['main']
                mainAnnoRow = [geneId, geneName, geneType, geneClass, txId, txName, txType, mainFeature]
                if args.mode == 'RNA':
                    minorFeature = mainAnnoPeakDict[geneId]['minor']
                    detailFeature = ','.join(mainAnnoPeakDict[geneId]['detail'])
                    mainAnnoRow.append(minorFeature)
                    mainAnnoRow.append(detailFeature)
                else:
                    distance = mainAnnoPeakDict[geneId]['distance']
                    mainAnnoRow.append(str(distance))
                outputRow = peakRow + mainAnnoRow + extraAnnoRow
                outputRowList.append(outputRow)
        else:
            if args.mode == 'RNA':
                mainAnnoRow = ['na' for i in range(10)]
            else:
                mainAnnoRow = ['na' for i in range(9)]
            outputRow = peakRow + mainAnnoRow + extraAnnoRow
            outputRowList.append(outputRow)
    else:
        appendRow = ['intergenic', 'intergenic' 'intergenic', 'intergenic']
        appendRow.extend(['na' for i in range(11)])
        appendRow[7] = 'intergenic'
        outputRow = peakRow + appendRow
        outputRowList.append(outputRow)

# writing to output
with open(args.output, 'w') as out:
    for outputRow in outputRowList:
        line = '\t'.join(outputRow) + '\n'
        out.write(line)
