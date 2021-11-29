#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
import numpy as np
## own module
import bedutils

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-b', '--bed', action='store', type=str, required=True,
                    help='gene annotation file in bed12 format (4th column, gene_id:gene_name:gene_type:tx_id:tx_name:tx_type)')
parser.add_argument('-f', '--fdr', action='store', type=float,
                    default=0.1,
                    help='the cutoff for FDR ( < --fdr)')
parser.add_argument('-i', '--input', action='store', type=str, required=True,
                    help='the input JC|JCEC results file from rmats.py')
parser.add_argument('-l', '--lncdiff', action='store', type=float,
                    default=0.1,
                    help='the cutoff for IncLevelDifference ( | IncLevelDifference | >= --lncdiff)')
parser.add_argument('-p', '--pval', action='store', type=float,
                    default=0.05,
                    help='the cutoff for PValue ( < --pval)')
parser.add_argument('-e', '--exon', action='store', type=int, choices=[2, 3],
                    default=3,
                    help='number of exon types (short,long|short,median,long)')
parser.add_argument('-t', '--type', action='store', type=str, choices=['A5SS', 'A3SS', 'SE', 'MXE', 'RI'],
                    default='SE',
                    help='type of alternative splicing')
parser.add_argument('-o', '--output', action='store', type=str, required=True,
                    help='output result')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

geneDict = defaultdict(dict)
txDict = defaultdict(dict)
exonDict = defaultdict(dict)
with open(args.bed, 'r') as f:
    for line in f:
        if bool(re.match(r'^#', line)):
            continue
        row = line.strip().split('\t')
        bedobj = bedutils.buildbed(row)
        tx2feature = bedobj.decode()
        ## get chroms, strand
        chrom = tx2feature.chr
        strand = tx2feature.strand
        ## get geneid
        infoList = tx2feature.name.split(':')
        geneid = infoList[0]
        geneName  = infoList[1]
        geneType = infoList[2]
        txid = infoList[3]
        txName = infoList[4]
        txType = infoList[5]
        ## store gene information
        if geneid not in geneDict:
            geneDict[geneid]['chrom'] = chrom
            geneDict[geneid]['strand'] = strand
            geneDict[geneid]['geneName'] = geneName
            geneDict[geneid]['geneType'] = geneType
            geneDict[geneid]['txid'] = defaultdict(dict)
            geneDict[geneid]['exon'] = defaultdict(dict)
        geneDict[geneid]['txid'][txid]['feature'] = tx2feature
        geneDict[geneid]['txid'][txid]['txName'] = txName
        geneDict[geneid]['txid'][txid]['txType'] = txType
        ## get exons
        exons = tx2feature.exon
        for exon in exons:
            exonCor = '_'.join(map(str, exon))
            if exonCor not in geneDict[geneid]['exon']:
                geneDict[geneid]['exon'][exonCor] = defaultdict(dict)
            geneDict[geneid]['exon'][exonCor]['txid'][txid] = ''

## calculate the median length of exons of all genes
exonLenList = []
for geneid in sorted(geneDict.keys()):
    for exonCor in sorted(geneDict[geneid]['exon'].keys()):
        exon = list(map(int, exonCor.split('_')))
        exonLen = exon[1] - exon[0]
        exonLenList.append(exonLen)

exonLenQ1 = np.percentile(exonLenList, 25)
exonLenQ3 = np.percentile(exonLenList, 75)
exonLenMedian = np.median(exonLenList)

for geneid in sorted(geneDict.keys()):
    exonCors = sorted(geneDict[geneid]['exon'].keys())
    ## decode the exon list and sorted by ascending exon start and then exon end
    #exonList = sorted(map(lambda x: list(map(int, x.split('_'))), exonCors),  key = lambda x:[x[0], x[1]])
    ### calculate the length of each exons
    #exonLens = list(map(lambda x: x[1] - x[0], exonList))
    ### determin the median of exon length
    #exonLenQ1 = np.percentile(exonLens, 25)
    #exonLenQ3 = np.percentile(exonLens, 75)
    #geneDict[geneid]['exonLenQ1'] = exonLenQ1
    #geneDict[geneid]['exonLenQ3'] = exonLenQ3
    ## determin the position classification of exons in each isoform
    ## determin the exon length classification, short ( length < exonLenQ1)
    ## median ( exonLenQ1 <= length <= exonLenQ3), long (length > exonLenQ3)
    for exonCor in exonCors:
        ## determin the exon length classification
        exon = list(map(int, exonCor.split('_')))
        exonLen = exon[1] - exon[0]
        if args.exon == 2:
            if exonLen < exonLenMedian:
                lenClass = 'short'
            elif exonLen >= exonLenMedian:
                lenClass = 'long'
        else:
            if exonLen < exonLenQ1:
                lenClass = 'short'
            elif exonLen > exonLenQ3:
                lenClass = 'long'
            else:
                lenClass = 'median'
        geneDict[geneid]['exon'][exonCor]['lenClass'] = lenClass
        ##tag the position class of exon in each isoforms
        for txid in sorted(geneDict[geneid]['exon'][exonCor]['txid'].keys()):
            tx2feature = geneDict[geneid]['txid'][txid]['feature']
            txidExons = tx2feature.exon
            txidExonNum = len(txidExons)
            for i in range(txidExonNum):
                if exon == txidExons[i]:
                    if i == 0:
                        geneDict[geneid]['exon'][exonCor]['txid'][txid] = 'first'
                        break
                    elif i == (txidExonNum - 1):
                        geneDict[geneid]['exon'][exonCor]['txid'][txid] = 'last'
                        break
                    else:
                        geneDict[geneid]['exon'][exonCor]['txid'][txid] = 'internal'
                        break

with open(args.input, 'r') as f, open(args.output, 'w') as out:
    row = f.readline().strip().split('\t')
    if args.type == 'SE':
        appendHeaderList = ['geneType', 'IJC_SAMPLE_1_Average', 'SJC_SAMPLE_1_Average', 'IJC_SAMPLE_2_Average', 'SJC_SAMPLE_2_Average', 'IncLevel1Ave', 'IncLevel2Ave']
        appendHeaderList += ['skip_exon_len', 'upstream_exon_len', 'downstream_exon_len', 'inclusion_exons_length', 'skipped_exons_length'] 
        appendHeaderList += ['skip_exon_length_class', 'upstream_exon_length_class', 'downstream_exon_length_class'] 
        appendHeaderList += ['inclusion_isoforms', 'skip_isoforms', 'upstream_exon_position_in_inclusion_isoforms', 'down_exon_position_in_inclusion_isoforms\n']
    row = row + appendHeaderList
    out.write('\t'.join(row))
    for line in f:
        row = line.strip().split('\t')
        row[1] = row[1].replace('"', '')
        geneid = row[1]
        row[2] = row[2].replace('"', '')
        geneType = geneDict[geneid]['geneType']
        strand = geneDict[geneid]['strand']
        if args.type == 'SE':
            ## filter data with cutoff
            pval = float(row[18])
            fdr = float(row[19])
            lncLevelDiff = float(row[-1])
            if pval < args.pval and fdr < args.fdr and abs(lncLevelDiff) >= args.lncdiff:
                pass
            else:
                continue
            if bool(re.search(r'NA', row[20])):
                continue
            if bool(re.search(r'NA', row[21])):
                continue
            #IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2
            incJcS1Ave = np.mean(list(map(int, row[12].split(','))))
            sJcS1Ave = np.mean(list(map(int, row[13].split(','))))
            incJcS2Ave = np.mean(list(map(int, row[14].split(','))))
            sJcS2Ave = np.mean(list(map(int, row[15].split(','))))
            incLevel1Ave = np.mean(list(map(float, row[20].split(','))))
            incLevel2Ave = np.mean(list(map(float, row[21].split(','))))
            ## exon Upstream
            exonUpStart = row[7]
            exonUpEnd = row[8]
            exonUpCor = '_'.join([exonUpStart,exonUpEnd])
            ## exon Skip
            exonSkipStart = row[5]
            exonSkipEnd = row[6]
            exonSkipCor = '_'.join([exonSkipStart,exonSkipEnd])
            ## exon Downstream
            exonDownStart = row[9]
            exonDownEnd = row[10]
            exonDownCor = '_'.join([exonDownStart,exonDownEnd])
            ## exon length
            exonUpLen = int(row[8]) - int(row[7])
            exonSkipLen = int(row[6]) - int(row[5])
            exonDownLen = int(row[10]) - int(row[9])
            incExonTotalLength = exonUpLen + exonSkipLen + exonDownLen
            excExonTotalLength = exonUpLen + exonDownLen
            ## isoforms contain these exons
            exonUpTxidSet = set(geneDict[geneid]['exon'][exonUpCor]['txid'].keys())
            exonSkipTxidSet = set(geneDict[geneid]['exon'][exonSkipCor]['txid'].keys())
            exonDownTxidSet = set(geneDict[geneid]['exon'][exonDownCor]['txid'].keys())
            ## exon length class
            exonUpLenClass = geneDict[geneid]['exon'][exonUpCor]['lenClass']
            exonSkipLenClass = geneDict[geneid]['exon'][exonSkipCor]['lenClass']
            exonDownLenClass = geneDict[geneid]['exon'][exonDownCor]['lenClass']
            ## exon inclusion isoforms
            incTxidList = sorted(exonUpTxidSet.intersection(exonSkipTxidSet).intersection(exonDownTxidSet))
            incTxids = ','.join(incTxidList)
            if incTxids == "":
                incTxids = "unkown"
            ## get the position class of the upstream or downstream exon in inclusion isoforms
            upExonPosClassList = list(map(lambda x:geneDict[geneid]['exon'][exonUpCor]['txid'][x], incTxidList))
            upExonPosClasses = ','.join(upExonPosClassList)
            if upExonPosClasses == "":
                upExonPosClasses = "NA"
            downExonPosClassList = list(map(lambda x:geneDict[geneid]['exon'][exonDownCor]['txid'][x], incTxidList))
            downExonPosClasses = ','.join(downExonPosClassList)
            if downExonPosClasses == "":
                downExonPosClasses = "NA"
            ## exon skip isoforms
            excTxidList = sorted(exonUpTxidSet.intersection(exonDownTxidSet).difference(exonSkipTxidSet))
            excTxids = ','.join(excTxidList)
            if excTxids == "":
                excTxids = "unkown"
            ## output row
            if strand == '-':
                row += [geneType, incJcS1Ave, sJcS1Ave, incJcS2Ave, sJcS2Ave, incLevel1Ave, incLevel2Ave]
                row += [exonSkipLen, exonDownLen, exonUpLen, incExonTotalLength]
                row += [excExonTotalLength, exonSkipLenClass, exonDownLenClass]
                row += [exonUpLenClass, incTxids, excTxids, downExonPosClasses, upExonPosClasses + '\n']
            else:
                row += [geneType, incJcS1Ave, sJcS1Ave, incJcS2Ave, sJcS2Ave, incLevel1Ave, incLevel2Ave]
                row += [exonSkipLen, exonUpLen, exonDownLen, incExonTotalLength]
                row += [excExonTotalLength, exonSkipLenClass, exonUpLenClass]
                row += [exonDownLenClass, incTxids, excTxids, upExonPosClasses, downExonPosClasses + '\n']
        out.write('\t'.join(map(str, row)))
