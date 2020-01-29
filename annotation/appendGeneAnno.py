#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
import tempfile
import copy

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-anno', action='store', type=str, required=True,
                    help='Gene annotation file in bed12 format (main annotation) \
                    (4th column [gene_id:gene_name:gene_type:tx_id:tx_name:tx_type])')
parser.add_argument('-geneClassFile', action='store', type=str,
                    help='geneType-geneClass pairwise file')
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input peak file (bed6 or bed6+)')
parser.add_argument('-idType', action='store', type=str, choices=['gene', 'tx'],
                    default='gene',
                    help='input id type')
parser.add_argument('-idCol', action='store', type=int, required=True,
                    help='0-based index of id column')
parser.add_argument('-inCol', action='store', type=int, required=True,
                    help='0-based index of column of appended information (-1 means last column)')
parser.add_argument('-skip', action='store', type=int,
                    default=0,
                    help='skip # first line')
parser.add_argument('-noheader', action='store_true',
                    default=False,
                    help='no header line insert')
parser.add_argument('-onlyName', action='store_true',
                    default=False,
                    help='append gene_name only')
parser.add_argument('-nover', action='store_true',
                    default=False,
                    help='remove gene version')
parser.add_argument('-ncbiGeneInfo', action='store', type=str,
                    help='*.gene_info file downloaded from NCBI')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output result (coordiantes in 0-base)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

# built up gene-type pairwise relationships
if args.onlyName is False:
    geneClassDict = defaultdict(dict)
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

# built up tx-infor pairwise relationships
idDict = defaultdict(dict)
with open(args.anno) as f:
    for line in f:
        row = line.strip().split('\t')
        txInfo = row[3].split(':')
        if args.idType == 'gene':
            keyId = txInfo[0]
        else:
            keyId = txInfo[3]
        if args.nover:
            keyId = keyId.split('.')[0]
        idDict[keyId]['geneId'] = txInfo[0]
        idDict[keyId]['geneName'] = txInfo[1]
        idDict[keyId]['geneType'] = txInfo[2]
        idDict[keyId]['txName'] = txInfo[4]
        idDict[keyId]['txType'] = txInfo[5]

if args.noheader is False:
    if args.idType == 'gene':
        headerRow = ['GeneName', 'Synonyms', 'Description', 'GeneType', 'GeneClass']
    else:
        headerRow = ['GeneId', 'Synonyms', 'Description', 'GeneName', 'GeneType', 'GeneClass', 'TxName', 'TxType']
    if args.onlyName:
        headerRow = ['GeneName']
count = 1
with open(args.input, 'r') as f, open(args.output, 'w') as out:
    __ = [f.readline() for x in range(args.skip)]
    for line in f:
        row = line.rstrip().split('\t')
        if args.noheader is False and count == 1:
            extRow = headerRow
        else:
            keyId = row[args.idCol]
            if args.nover:
                keyId = keyId.split('.')[0]
            if keyId not in idDict:
                idDict[keyId]['geneId'] = 'na'
                idDict[keyId]['geneName'] = 'na'
                idDict[keyId]['geneType'] = 'na'
                idDict[keyId]['txName'] = 'na'
                idDict[keyId]['txType'] = 'na'
            if args.onlyName:
                geneName = idDict[keyId]['geneName']
                extRow = [geneName]
            else:
                geneName = idDict[keyId]['geneName']
                geneType = idDict[keyId]['geneType']
                geneClass = geneClassDict[geneType] if (geneType in geneClassDict) else 'Unkown'
                if args.idType == 'gene':
                    geneId = keyId
                else:
                    geneId = idDict[keyId]['geneId']
                keyIdNoVer = geneId.split('.')[0]
                if keyIdNoVer in ncbiGeneInfoDict:
                    synonyms, description = ncbiGeneInfoDict[keyIdNoVer]
                else:
                    synonyms, description = ['na', 'na']
                if args.idType == 'gene':
                    extRow = [geneName, synonyms, description, geneType, geneClass]
                else:
                    txName = idDict[keyId]['txName']
                    txType = idDict[keyId]['txType']
                    extRow = [geneId, geneName, synonyms, description, geneType, geneClass, txName, txType]
        if len(row) <= (args.inCol + 1) or args.inCol == -1:
            row.extend(extRow)
        else:
            row = row[0:args.inCol+1] + extRow + row[args.inCol+1:]
        out.write('\t'.join(row) + '\n')
        count += 1
