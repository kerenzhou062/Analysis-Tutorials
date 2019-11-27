#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
import gffutils
from PubAlbum import Anno, GeneClass

#usage: isoformAnno.py or isoformAnno.py -inpupt test.bed12 -format gtf --longest -rule exon -sqlite hg38v32 -type protein_coding -output test.bed12

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-extend', nargs='*', type=str,
                    help='extended gene types')
parser.add_argument('-format', action='store', type=str,
                    default='gtf', help='sqlite-file source: GFF3|gtf')
parser.add_argument('-input', action='store', type=str,
                     help='input bed12 file')
parser.add_argument('--longest', action='store_true',
                    default=False, help='turn on longest isoform mode')
parser.add_argument('-output', action='store', type=str,
                    default='test.bed12', help='output bed12')
parser.add_argument('-rule', action='store', type=str,
                    default='exon', help='rule for longest transcript isoform (exon|full)')
parser.add_argument('-sqlite', action='store', type=str,
                    default='hg38v32', help='sqlite3-file based databse from GFF3|gtf (eg. hg38v32)')
parser.add_argument('-type', action='store', type=str,
                    default='all', help='gene types for filtering \
                    (all|protein_coding|lncRNA|sncRNA|pseudogene)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

##public arguments
geneClass = GeneClass().filter(args.type, args.extend)
anno = Anno()
sqlite = anno.sqlite(args.sqlite, args.format)

geneDict = defaultdict(dict)
txDict = defaultdict(str)
sqliteDb = gffutils.FeatureDB(sqlite, keep_order=True)
with open(args.input, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        txId = row[3]
        txChildren = sqliteDb.children(id=txId)
        tx = next(txChildren)
        geneId = tx.attributes['gene_id'][0]
        geneType = tx.attributes['gene_type'][0]
        txType = tx.attributes['transcript_type'][0]
        #filter gene with specific gene class
        if args.type == 'all':
            pass
        elif geneType in geneClass and txType in geneClass:
            pass
        else:
            continue
        if args.longest:
            txLength = 0
            if args.rule == 'full':
                txLength = tx.end - tx.start + 1
                txChildren.close()
            while True:
                try:
                    children = next(txChildren)
                    if children.featuretype == 'exon':
                        txLength += children.end - children.start + 1
                except StopIteration:
                    break
            if geneId in geneDict:
                if txLength > geneDict[geneId]['len']:
                    preTxId = geneDict[geneId]['tx']
                    txDict[txId] = 1
                    txDict[preTxId] = 0
                    geneDict[geneId]['tx'] = txId
                    geneDict[geneId]['len'] = txLength
                else:
                    txDict[txId] = 0
            else:
                geneDict[geneId]['tx'] = txId
                geneDict[geneId]['len'] = txLength
                txDict[txId] = 1
        else:
            txDict[txId] = 1

with open(args.input, 'r') as f, open(args.output, 'w') as out:
    for line in f:
        row = line.split('\t')
        txId = row[3]
        if txId in txDict:
            if txDict[txId] == 0:
                continue
            tx = sqliteDb[txId]
            geneId = tx.attributes['gene_id'][0]
            geneName = tx.attributes['gene_name'][0]
            geneType = tx.attributes['gene_type'][0]
            txName = tx.attributes['transcript_name'][0]
            txType = tx.attributes['transcript_type'][0]
            row[3] = ':'.join([geneId, geneName, geneType, 
                txId, txName, txType])
            out.write('\t'.join(row))
