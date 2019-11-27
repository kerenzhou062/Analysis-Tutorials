#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
import datetime
from multiBioPro import BedMan

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                    help='The integrated mirTarget-seq file')
parser.add_argument('-cds', action='store_true',
                    default=False, help='cds flag')
parser.add_argument('-circRNA', action='store_true',
                    default=False, help='cds flag')
parser.add_argument('-type', action='store', type=str,
                    default='bed12', help='bed6 or bed12(default:bed12)')
parser.add_argument('-output', action='store', type=str,
                    default="RBPanno.txt", help='The output clusterID file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

starttime = datetime.datetime.now()

RBPgeneExpDict = defaultdict(dict)
clusterDict = defaultdict(dict)
if args.type == 'bed12':
    with open(args.input, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            clusterID = row[3]
            start = int(row[1])
            end = int(row[2])
            bed12Row = row[11:]
            # [[exonblock], [intronblock], [thickup, thick, thickdown]]
            # or [[exonblock], [intronblock]]
            if start < int(bed12Row[1]) or end > int(bed12Row[2]):
                continue
            decodeList = BedMan.decodeBed12(bed12Row)
            locus = [start, end]
            if len(decodeList) == 2:
                if args.circRNA:
                    pass
                elif args.cds:
                    continue
                else:
                    pass
            else:
                if args.cds:
                    pass
                else:
                    continue
            txInfo = bed12Row[3].split('|')
            txID = txInfo[0]
            txName = txInfo[1]
            geneID = txInfo[3]
            geneName = txInfo[4]
            geneType = txInfo[5]
            strand = bed12Row[5]
            if strand == '+':
                thickTypeList = ['5\'UTR', 'CDS', '3\'UTR']
            else:
                thickTypeList = ['3\'UTR', 'CDS', '5\'UTR']
            exonBlock = decodeList[0]
            intronBlock = decodeList[1]
            exonBlockOverlap = list(map(lambda x:BedMan.overlap(locus, x),
                exonBlock))
            exonOverlapList = list()
            exonCount = len(exonBlockOverlap)
            for i in range(exonCount):
                if exonBlockOverlap[i]:
                    if strand == '+':
                        exonNum = 'Exon-' + str(i + 1)
                    else:
                        exonNum = 'Exon-' + str(exonCount - i)
                    exonOverlapList.append(exonNum)
            exonBlockLocate = ','.join(exonOverlapList)
            # design for cds
            if args.cds:
                if len(decodeList) == 3:
                    thickBlock = decodeList[2]
                    thickBlockOverlapList = list()
                    for i in thickBlock:
                        if i:
                            overlap = list(map(lambda x:BedMan.overlap(locus, x),i))
                            thickBlockOverlapList.append(overlap)
                        else:
                            thickBlockOverlapList.append([])
                    thickOverlapList = list()
                    for i in range(3):
                        if sum(thickBlockOverlapList[i]):
                            thickOverlapList.append(thickTypeList[i])
                    thickBlockLocate = ','.join(thickOverlapList)
                else:
                    thickBlockLocate = 'Exon'
            else:
                thickBlockLocate = 'Exon'
            if clusterID not in clusterDict:
                clusterDict[clusterID]['info'] = row[0:11]
                clusterDict[clusterID]['gene'] = defaultdict(list)
                clusterDict[clusterID]['gene'][geneID].append([geneName, geneType, txID,
                    txName, exonBlockLocate, thickBlockLocate])
            else:
                clusterDict[clusterID]['gene'][geneID].append([geneName, geneType, txID,
                    txName, exonBlockLocate, thickBlockLocate])
            # RBP-gene statistics
            RBP = row[8]
            clipIDList = row[10].split(',')
            clipExpList = list(map(lambda x: x.split('-')[0], clipIDList))
            if RBP not in RBPgeneExpDict:
                RBPgeneExpDict[RBP][geneID] = clipExpList
            else:
                if geneID not in RBPgeneExpDict[RBP]:
                    RBPgeneExpDict[RBP][geneID] = clipExpList
                else:
                    RBPgeneExpDict[RBP][geneID].extend(clipExpList)
else:
    with open(args.input, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            clusterID = row[3]
            start = int(row[1])
            end = int(row[2])
            bed6Row = row[11:]
            if start < int(bed6Row[1]) or end > int(bed6Row[2]):
                continue
            geneInfo = bed6Row[3].split('|')
            geneID = geneInfo[3]
            geneName = geneInfo[4]
            geneType = geneInfo[5]
            if clusterID not in clusterDict:
                clusterDict[clusterID]['info'] = row[0:11]
                clusterDict[clusterID]['gene'] = list()
                clusterDict[clusterID]['gene'].append([geneID, geneName, geneType])
            else:
                clusterDict[clusterID]['gene'].append([geneID, geneName, geneType])
            # RBP-gene statistics
            RBP = row[8]
            clipIDList = row[10].split(',')
            clipExpList = list(map(lambda x: x.split('-')[0], clipIDList))
            if RBP not in RBPgeneExpDict:
                RBPgeneExpDict[RBP][geneID] = clipExpList
            else:
                if geneID not in RBPgeneExpDict[RBP]:
                    RBPgeneExpDict[RBP][geneID] = clipExpList
                else:
                    RBPgeneExpDict[RBP][geneID].extend(clipExpList)

for RBP in sorted(RBPgeneExpDict.keys()):
    for geneID in sorted(RBPgeneExpDict[RBP].keys()):
        RBPgeneExpDict[RBP][geneID] = str(len(set(RBPgeneExpDict[RBP][geneID])))

if args.type == 'bed12':
    if args.circRNA:
        headerList = ['lineID', 'chromosome', 'narrowStart', 'narrowEnd', 'clusterID',
            'clusterClipExpNum', 'strand', 'broadStart', 'broadEnd',
            'RBP', 'clipExpNum', 'clusterClipSiteNum', 'clipID', 'geneID',
            'geneName', 'geneType', 'txIDcat', 'txInfo']
    else:
        headerList = ['lineID', 'chromosome', 'narrowStart', 'narrowEnd', 'clusterID',
            'clusterClipExpNum', 'strand', 'broadStart', 'broadEnd',
            'RBP', 'clipExpNum', 'clusterClipSiteNum', 'clipID', 'geneID',
            'geneName', 'geneType', 'txInfo']
else:
    headerList = ['lineID', 'chromosome', 'narrowStart', 'narrowEnd', 'clusterID',
        'clusterClipExpNum', 'strand', 'broadStart', 'broadEnd',
        'RBP', 'clipExpNum', 'clusterClipSiteNum', 'clipID', 'geneID',
        'geneName', 'geneType']

if args.type == 'bed12':
    with open(args.output, 'w') as out:
        out.write('\t'.join(headerList) + '\n')
        clusterIDs = sorted(clusterDict.keys())
        count = 1
        for clusterID in clusterIDs:
            tempDict = clusterDict[clusterID]
            infoList = tempDict['info']
            RBP = infoList[8]
            geneIDlist = sorted(tempDict['gene'].keys())
            for geneID in geneIDlist:
                # add clipExpNum and degraExpNum to clusterInfo
                tempInfoList = infoList[0:]
                clipExpNum = RBPgeneExpDict[RBP][geneID]
                tempInfoList[9] = '\t'.join([clipExpNum, tempInfoList[9]])
                tempInfo = '\t'.join(tempInfoList)
                # cat transcript info
                tempTxList = tempDict['gene'][geneID]
                geneName = tempTxList[0][0]
                geneType = tempTxList[0][1]
                tempList = [geneID, geneName, geneType]
                tempTxCatList = list()
                tempTxInfoCatList = list()
                for tx in tempTxList:
                    txInfo = ':'.join(tx[2:])
                    tempTxInfoCatList.append(txInfo)
                    txID = tx[2]
                    tempTxCatList.append(txID)
                if args.circRNA:
                    tempTxCat = ','.join(sorted(set(tempTxCatList)))
                    tempList.append(tempTxCat)
                tempTx = '|'.join(tempTxInfoCatList)
                tempList.append(tempTx)
                temp = '\t'.join(tempList) + '\n'
                out.write('\t'.join([str(count), tempInfo, temp]))
                count += 1
else:
    with open(args.output, 'w') as out:
        out.write('\t'.join(headerList) + '\n')
        clusterIDs = sorted(clusterDict.keys())
        count = 1
        for clusterID in clusterIDs:
            tempDict = clusterDict[clusterID]
            infoList = tempDict['info']
            RBP = infoList[8]
            for geneInfo in sorted(tempDict['gene'], key=lambda x:x[0]):
                geneID = geneInfo[0]
                geneName = geneInfo[1]
                geneType = geneInfo[2]
                # add clipExpNum and degraExpNum to clusterInfo
                tempInfoList = infoList[0:]
                clipExpNum = RBPgeneExpDict[RBP][geneID]
                tempInfoList[9] = '\t'.join([clipExpNum, tempInfoList[9]])
                tempInfo = '\t'.join(tempInfoList)
                out.write('\t'.join([str(count), tempInfo, geneID, geneName, geneType]) + '\n')
                count += 1

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
