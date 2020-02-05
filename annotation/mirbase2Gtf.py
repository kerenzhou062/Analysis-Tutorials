#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict

#usage: mirBase2Bed.py -inpupt miRNA.v22.gff3 --tx -output miRNA.v22.bed12

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str,
                     help='input GFF3 file')
parser.add_argument('-mode', action='store', type=str,
                    help='primary|mature mode')
parser.add_argument('-output', action='store', type=str,
                    help='output gtf')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

if args.mode not in ['primary', 'mature']:
    sys.stderr('Wrong -mode!')
    parser.print_help()
    parser.exit()

# primary mode: gene:primary, transcript:primary_1, exon:primary_1_1
# mature mode: gene:primary, transcript:miRNA, exon:miRNA_1

geneType = 'miRNA'
attDict = defaultdict(dict)
with open(args.input, 'r') as f:
    for line in f:
        if re.match(r'^#', line):
            continue
        row = line.strip().split('\t')
        feature = row[2]
        attTempDict = defaultdict(dict)
        for item in row[8].split(';'):
            key = item.split('=')[0]
            value = item.split('=')[1]
            attTempDict[key] = value
        attDict[attTempDict['ID']] = attTempDict

with open(args.input, 'r') as f, open(args.output, 'w') as out:
    for line in f:
        if re.match(r'^#', line):
            continue
        row = line.strip().split('\t')
        feature = row[2]
        chrom = row[0]
        start = row[3]
        end = row[4]
        strand = row[6]
        attTempDict = defaultdict(dict)
        for item in row[8].split(';'):
            key = item.split('=')[0]
            value = item.split('=')[1]
            attTempDict[key] = value
        # output
        if feature == 'miRNA_primary_transcript':
            row[2] = 'gene'
            geneId = attTempDict['ID']
            geneName = attTempDict['Name']
            row[-1] = 'gene_id "{geneId}"; gene_type "{geneType}"; gene_name "{geneName}";'
            row[-1] = row[-1].format(**vars())
            out.write('\t'.join(row) + '\n')
            if args.mode == 'primary':
                # gene:primary, transcript:primary_1, exon:primary_1_1
                row[2] = 'transcript'
                row[-1] += ' transcript_id "{geneId}_1"; '
                row[-1] += 'transcript_type "{geneType}"; transcript_name "{geneName}_1";'
                row[-1] = row[-1].format(**vars())
                out.write('\t'.join(row) + '\n')
                row[2] = 'exon'
                row[-1] += ' exon_number "1"; exon_id "{geneId}_1_1";'
                row[-1] = row[-1].format(**vars())
                out.write('\t'.join(row) + '\n')
        elif feature == 'miRNA':
            mirnaId = attTempDict['ID']
            mirnaName = attTempDict['Name']
            parentId = attTempDict['Derives_from']
            geneId = parentId
            geneName = attDict[geneId]['Name']
            row[-1] = 'gene_id "{geneId}"; gene_type "{geneType}"; gene_name "{geneName}";'
            if args.mode == 'mature':
                #
                row[2] = 'transcript'
                row[-1] += ' transcript_id "{mirnaId}";'
                row[-1] += ' transcript_type "{geneType}"; transcript_name "{mirnaName}";'
                row[-1] = row[-1].format(**vars())
                out.write('\t'.join(row) + '\n')
                row[2] = 'exon'
                row[-1] += ' transcript_id "{mirnaId}";'
                row[-1] += ' transcript_type "{geneType}"; transcript_name "{mirnaName}";'
                row[-1] += ' exon_number "1"; exon_id "{mirnaId}_1";'
                row[-1] = row[-1].format(**vars())
                out.write('\t'.join(row) + '\n')
