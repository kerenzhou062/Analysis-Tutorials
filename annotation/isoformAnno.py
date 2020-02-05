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
                    default='exon', help='rule for longest transcript isoform (exon|full), work with --longest')
parser.add_argument('-sqlite', action='store', type=str,
                    help='sqlite3-file based databse from GFF3|gtf (eg. hg38v32)')
parser.add_argument('-gtf', action='store', type=str,
                    help='gtf annotation file, not work with --longest')
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
geneDict = defaultdict(dict)
txDict = defaultdict(dict)

if bool(args.sqlite):
    sqlite = anno.sqlite(args.sqlite, args.format)
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
else:
    with open(args.gtf, 'r') as f:
        for line in f:
            if bool(re.match(r'^#', line)):
                continue
            row = line.strip().split('\t')
            fType = row[2]
            if fType != 'exon':
                continue
            txId = re.findall(r'transcript_id "(.+?)";', row[-1])[0]
            if txId not in txDict:
                geneId = re.findall(r'gene_id "(.+?)";', row[-1])[0]
                geneName = re.findall(r'gene_name "(.+?)";', row[-1])[0]
                geneType = re.findall(r'gene_type "(.+?)";', row[-1])[0]
                txName = re.findall(r'transcript_name "(.+?)";', row[-1])[0]
                txType = re.findall(r'transcript_type "(.+?)";', row[-1])[0]
                txDict[txId]['gene_id'] = geneId
                txDict[txId]['gene_name'] = geneName
                txDict[txId]['gene_type'] = geneType
                txDict[txId]['transcript_name'] = txName
                txDict[txId]['transcript_type'] = txType

with open(args.input, 'r') as f, open(args.output, 'w') as out:
    for line in f:
        row = line.split('\t')
        txId = row[3]
        if bool(args.sqlite):
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
        else:
            if txId in txDict:
                geneId = txDict[txId]['gene_id']
                geneName = txDict[txId]['gene_name']
                geneType = txDict[txId]['gene_type']
                txName = txDict[txId]['transcript_name']
                txType = txDict[txId]['transcript_type']
                row[3] = ':'.join([geneId, geneName, geneType, 
                    txId, txName, txType])
                out.write('\t'.join(row))
