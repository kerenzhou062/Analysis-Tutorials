#!/usr/bin/env python3
import os
import sys
import argparse
import tempfile
import subprocess
import re
from collections import defaultdict
# in-house module
import bedutils

parser = argparse.ArgumentParser(
    description="This script is used for annotating bed output from convergeAaToSegBed.py",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', action='store', type=str,
                    required=True,
                    help='bed output from convergeAaToSegBed.py')
parser.add_argument('-g', '--geneClassFile', action='store', type=str,
                    help='geneType-geneClass pairwise file')
parser.add_argument('--ncbiGeneInfo', action='store', type=str,
                    help='*.gene_info file downloaded from NCBI')
parser.add_argument('-o', '--output', action='store', type=str,
                    required=True,
                    help='output bed')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--bed12', nargs='+', type=str,
                    help='Gene annotation file in bed12 format \
                    (4th column [gene_id:gene_name:gene_type:tx_id:tx_name:tx_type])')
group.add_argument('--bed6', nargs='+', type=str,
                    help='Gene annotation file in bed6 format  \
                    (4th column[gene_name:gene_type:attribute:gene_class])')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

# function
def decodeInclude(locusA, locusB, atype):
    ctypeDict = defaultdict(dict)
    ctypeDict['.'] = {0:'whole', 1:'right', -1:'left', 2:'overlay'}
    ctypeDict['+'] = {0:'whole', 1:'TTS', -1:'TSS', 2:'overlay'}
    ctypeDict['-'] = {0:'whole', 1:'TSS', -1:'TTS', 2:'overlay'}
    ## use bedutils
    bedops = bedutils.bedops(locusA, locusB)
    include = bedops.include()
    ctype = include.ctype
    uniqSegId = bedops.a.name
    segStrand = bedops.a.strand
    annoStrand = bedops.b.strand
    annoLength = bedops.b.length
    otype = ctypeDict[annoStrand][ctype]
    if annoStrand == '.':
        orient = '.'
    elif segStrand == annoStrand:
        orient = '+'
    else:
        orient = '-'
    cloverh = bedops.cloverh
    croverh = bedops.croverh
    absoverh = abs(cloverh) + abs(croverh)
    overlap = bedops.intersect()
    fracA = overlap.ifracA
    fracB = overlap.ifracB
    overlapLength = overlap.ilength
    BLocus = [bedops.b.start, bedops.b.end, bedops.b.strand]
    includeList = BLocus +[ctype, otype, orient]
    includeList += [cloverh, croverh, absoverh, overlapLength, fracA, fracB]
    if atype == 'gene':
        annoInfoList = bedops.b.name.split('|')
        geneId = annoInfoList[0]
        if geneId not in includeDict[uniqSegId]:
            includeDict[uniqSegId][geneId] = defaultdict(dict)
        includeDict[uniqSegId][geneId]['gene'] = includeList
    else:
        annoInfoList = bedops.b.name.split(':')
        geneId = annoInfoList[0]
        txId = annoInfoList[3]
        if geneId not in includeDict[uniqSegId]:
            includeDict[uniqSegId][geneId] = defaultdict(dict)
        includeDict[uniqSegId][geneId]['tx'][txId] = includeList
    return 1

# main

# built up gene-type pairwise relationships
geneClassDict = defaultdict(dict)
if bool(args.geneClassFile):
    with open(args.geneClassFile) as f:
        for line in f:
            row = line.strip().split('\t')
            mainType = row[0]
            geneType = row[1]
            geneClassDict[geneType] = mainType

# built up gene-info pairwise relationships, if args.ncbiGeneInfo
ncbiGeneInfoDict = defaultdict(dict)
if bool(args.ncbiGeneInfo):
    with open(args.ncbiGeneInfo, 'r') as f:
        __ = f.readline()
        for line in f:
            row = line.strip().split('\t')
            dbXrefsDict = defaultdict(dict)
            ## MIM:138670|HGNC:HGNC:5|Ensembl:ENSG00000121410
            for dbXrefs in row[5].split('|'):
                db = dbXrefs.split(':')[0]
                geneId = dbXrefs.split(':')[-1]
                dbXrefsDict[db] = geneId
            if 'Ensembl' in dbXrefsDict:
                geneId = dbXrefsDict['Ensembl']
                synonyms = row[4] if row[4] != '-' else 'na'
                description = row[8] if row[8] != '-' else 'na'
                ncbiGeneInfoDict[geneId] = [synonyms, description]

ecdnaNameRow = list()
ecdnaDict = defaultdict(dict)
with open(args.input, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        if re.match(r'^#', line):
            ecdnaNameRow = row
        else:
            uniqSegId = row[3]
            ecdnaDict[uniqSegId] = row

if bool(args.bed12):
    annoBedList = args.bed12
else:
    annoBedList = args.bed6
annoInfoDict = defaultdict(dict)
annoInfoDict['tx'] = defaultdict(dict)
annoInfoDict['gene'] = defaultdict(dict)
geneBedDict = defaultdict(dict)
count = 1
for annoBed in annoBedList:
    with open(annoBed, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            strand = row[5]
            nameRow = row[3].split(':')
            if bool(args.bed12):
                geneId, geneName, geneType, txId, txName, txType = nameRow
                if geneId not in geneBedDict:
                    info = '|'.join([geneId, geneName, geneType])
                    geneBedDict[geneId]['gene'] = [chrom, start, end, info, 1, strand]
                else:
                    geneBedDict[geneId]['gene'][4] += 1
                    if start < geneBedDict[geneId]['gene'][1]:
                        geneBedDict[geneId]['gene'][1] = start
                    if geneBedDict[geneId]['gene'][2] < end:
                        geneBedDict[geneId]['gene'][2] = end
                # build tx information dict
                annoInfoDict['tx'][txId]['geneId'] = geneId
                annoInfoDict['tx'][txId]['geneName'] = geneName
                annoInfoDict['tx'][txId]['geneType'] = geneType
                annoInfoDict['tx'][txId]['txType'] = txType
                annoInfoDict['gene'][geneId]['geneName'] = geneName
                annoInfoDict['gene'][geneId]['geneType'] = geneType
                if txType in geneClassDict:
                    annoInfoDict['tx'][txId]['geneClass'] = geneClassDict[txType]
                elif geneType in geneClassDict:
                    annoInfoDict['tx'][txId]['geneClass'] = geneClassDict[geneType]
                else:
                    annoInfoDict['tx'][txId]['geneClass'] = 'other'
    
                if geneType in geneClassDict:
                    annoInfoDict['gene'][geneId]['geneClass'] = geneClassDict[geneType]
                else:
                    annoInfoDict['gene'][geneId]['geneClass'] = 'other'
            else:
                geneName, geneType, attribute, geneClass = nameRow
                if re.search(r'miR-|let-', geneType, re.IGNORECASE):
                    geneId = '='.join([geneName, str(count)])
                    geneName = geneType
                    geneType = 'miRNA'
                    geneClass = 'miRNA'
                else:
                    if attribute != 'NA':
                        geneName = ':'.join([geneName, attribute])
                    geneId = '='.join([geneName, str(count)])
                row[3] = '|'.join([geneId, geneName, geneType])
                annoInfoDict['gene'][geneId]['geneName'] = geneName
                annoInfoDict['gene'][geneId]['geneType'] = geneType
                annoInfoDict['gene'][geneId]['geneClass'] = geneClass
                geneBedDict[geneId]['gene'] = row
            count += 1

annoBed6Tmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
with open(annoBed6Tmp.name, 'w') as temp:
    for geneId in sorted(geneBedDict.keys()):
        row = geneBedDict[geneId]['gene']
        bed6Line = '\t'.join(map(str, row)) + '\n'
        temp.write(bed6Line)

annoBed6Tmp.seek(0)


geneAnnoCommand = 'bedtools intersect -a {0} -b {1} -wa -wb'.format(annoBed6Tmp.name, args.input)
geneAnnoResList = bytes.decode(subprocess.check_output(geneAnnoCommand, shell=True)).split('\n')
annoBed6Tmp.close()

includeDict = defaultdict(dict)
for line in geneAnnoResList:
    if bool(line) is False:
        continue
    row = line.strip().split('\t')
    annoRow = row[0:6]
    ecdnaRow = row[6:]
    bedRow = ecdnaRow[0:6]
    decode = decodeInclude(bedRow, annoRow, 'gene')
if bool(args.bed12):
    txAnnoResList = list()
    for annoBed in annoBedList:
        txAnnoCommand = 'bedtools intersect -a {0} -b {1} -wa -wb'.format(annoBed, args.input)
        txAnnoResList += bytes.decode(subprocess.check_output(txAnnoCommand, shell=True)).split('\n')
    
    for line in txAnnoResList:
        if bool(line) is False:
            continue
        row = line.strip().split('\t')
        #txRow = row[0:12]
        annoRow = row[0:6]
        ecdnaRow = row[12:]
        bedRow = ecdnaRow[0:6]
        decode = decodeInclude(bedRow, annoRow, 'tx')

row = ecdnaNameRow + ["geneId", "geneName", "synonyms", "description", "geneClass", "geneType", "txId", "txType"]
row += ["annoStart", "annoEnd", "annoStrand", "ctype", "otype", "orient"] 
row += ["cloverh", "croverh", "absoverh", "overlapLength", "fracA", "fracB"]

overhangPri = ['whole', 'TTS', 'TSS', 'overlay']
with open(args.output, 'w') as out:
    out.write('\t'.join(row) + '\n')
    for uniqSegId in sorted(ecdnaDict.keys()):
        if uniqSegId in includeDict:
            geneIdList = sorted(includeDict[uniqSegId].keys())
            keptGeneDict = defaultdict(dict)
            for geneId in geneIdList:
                #[cloverh, croverh, absoverh, ctype, otype, orient]
                tempDict = includeDict[uniqSegId][geneId]
                geneInList = tempDict['gene']
                keptTxId = ''
                ## 'gene': [88781745, 88851628, '-', 0, 'whole', '-', 190324, -9386, 199710, 69883, 0.2592166710559992, 1.0]
                ## 'tx': {'ENST00000301015.14_4': [88781745, 88851628, '-', 0, 'whole', '-', 190324, -9386, 199710, 69883, 0.2592166710559992, 1.0]}
                if geneInList[4] != 'whole':
                    if 'tx' in tempDict:
                        for txId in sorted(tempDict['tx'].keys()):
                            if bool(keptTxId) is False:
                                keptTxId = txId
                            else:
                                if overhangPri.index(tempDict['tx'][txId][4]) < overhangPri.index(tempDict['tx'][keptTxId][4]):
                                    keptTxId = txId
                                else:
                                    if tempDict['tx'][txId][8] < tempDict['tx'][keptTxId][8]:
                                        keptTxId = txId
                if bool(keptTxId):
                    keptGeneDict[geneId]['include'] = tempDict['tx'][keptTxId]
                    keptGeneDict[geneId]['annoType'] = 'isoform'
                    keptGeneDict[geneId]['annoId'] = keptTxId
                else:
                    keptGeneDict[geneId]['include'] = geneInList
                    keptGeneDict[geneId]['annoType'] = 'gene'
                    keptGeneDict[geneId]['annoId'] = geneId
            geneIdList = sorted(keptGeneDict.keys())
            for geneId in geneIdList:
                include = keptGeneDict[geneId]['include']
                annoType = keptGeneDict[geneId]['annoType']
                annoId = keptGeneDict[geneId]['annoId']
                geneType = annoInfoDict['gene'][geneId]['geneType']
                geneName = annoInfoDict['gene'][geneId]['geneName']
                if annoType == 'gene':
                    txId = 'na'
                    txType = 'na'
                    geneClass = annoInfoDict['gene'][geneId]['geneClass']
                else:
                    txId = annoId
                    txType = annoInfoDict['tx'][txId]['txType']
                    geneClass = annoInfoDict['tx'][txId]['geneClass']
                ensemblId = geneId.split('.')[0]
                if ensemblId in ncbiGeneInfoDict:
                    synonyms = ncbiGeneInfoDict[ensemblId][0]
                    description = ncbiGeneInfoDict[ensemblId][1]
                else:
                    synonyms = 'na'
                    description = 'na'
                if bool(args.bed6):
                    realGeneId = geneId.split('=')[0]
                else:
                    realGeneId = geneId
                annoRow = [realGeneId, geneName, synonyms, description, geneClass, geneType, txId, txType]
                annoRow += include
                row = ecdnaDict[uniqSegId] + annoRow
                out.write('\t'.join(map(str, row)) + '\n')
        else:
            include = ['.', '.', '.', 'na', 'na', 'na', 0, 0, 0, 0, 0, 0, ]
            annoRow = ['intergenic', 'intergenic', 'na', 'na', 'intergenic', 'intergenic', 'na', 'na']
            annoRow += include
            row = ecdnaDict[uniqSegId] + annoRow
            out.write('\t'.join(map(str, row)) + '\n')
