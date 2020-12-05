#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-ctk', nargs='+', type=str,
                    help='final DIR (or file) of CTK results')
parser.add_argument('-rbsb', nargs='+', type=str,
                    help='final DIR (or file) of rbsSeeker(bowtie) results')
parser.add_argument('-rbss', nargs='+', type=str,
                    help='final DIR (or file) of rbsSeeker(STAR) results')
parser.add_argument('-grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('-grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('-prefix', action='store', type=str, required=True,
                    help='output prefix')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='output directory')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

# function
def ConsInputDict (aDict, aInput, key, regex, kept):
    if bool(aInput):
        finalList = list()
        inputList = aInput
        if os.path.isdir(inputList[0]):
            fileList = sorted(glob(os.path.join(inputList[0], '**', '*.bed'), recursive=True))
            for file in fileList:
                if os.path.isfile(file) is False:
                    continue
                fileName = os.path.split(file)[-1]
                if bool(regex):
                    if kept:
                        if bool(regex.search(fileName)) is False:
                            continue
                    else:
                        if bool(regex.search(fileName)) is True:
                            continue
                finalList.append(file)
        else:
            for file in fileList:
                if os.path.isfile(file):
                    finalList.append(file)
        aDict[key] = finalList

##public arguments
ctkBool = bool(args.ctk) is False
rbsbBool = bool(args.rbsb) is False
rbssBool = bool(args.rbss) is False
if ctkBool and rbsbBool and rbssBool:
    parser.print_help()
    parser.exit()

if bool(args.grepKept):
    kept = True
    regex = re.compile(r'{0}'.format(args.grepKept))
elif bool(args.grepExpel):
    kept = False
    regex = re.compile(r'{0}'.format(args.grepExpel))
else:
    regex = False

inputDict = defaultdict(list)
ConsInputDict(inputDict, args.ctk, 'CTK', regex, kept)
ConsInputDict(inputDict, args.rbsb, 'rbsSeeker_bowtie', regex, kept)
ConsInputDict(inputDict, args.rbss, 'rbsSeeker_STAR', regex, kept)

lRegex = re.compile(r'^#')
ctRegex = re.compile(r'rbsSeeker_CT|CIMS')
trRegex = re.compile(r'rbsSeeker_Truncation|CITS')

siteDict = defaultdict(dict)
keyList = sorted(inputDict.keys())
for key in keyList:
    for file in sorted(inputDict[key]):
        fileName = os.path.basename(file)
        rep = re.findall(r'rep\d+', fileName)[0]
        with open(file, 'r') as f:
            for line in f:
                if bool(lRegex.search(line)):
                    continue
                row = line.strip().split('\t')
                site = '\t'.join([row[0], row[1], row[2], row[5]])
                idCol = row[3]
                if bool(ctRegex.search(idCol)) and bool(trRegex.search(idCol)):
                    siteType = 'CIMS|CITS'
                elif bool(ctRegex.search(idCol)):
                    siteType = 'CIMS'
                elif bool(trRegex.search(idCol)):
                    siteType = 'CITS'
                else:
                    continue
                if key == 'CTK':
                    readCount = 0
                else:
                    readCount = int(row[4])
                seq = idCol.split('=')[0]
                if len(seq) != 5:
                    seq = 'none'
                if site in siteDict:
                    if siteDict[site]['read'] < readCount:
                        siteDict[site]['read'] = readCount
                    siteDict[site]['key'].append(key)
                    siteDict[site]['rep'].append(rep)
                    siteDict[site]['type'].append(siteType)
                    siteDict[site]['store'].append(idCol)
                    siteDict[site]['file'].append(fileName)
                else:
                    siteDict[site]['read'] = int(readCount)
                    siteDict[site]['key'] = [key]
                    siteDict[site]['rep'] = [rep]
                    siteDict[site]['type'] = [siteType]
                    siteDict[site]['store'] = [idCol]
                    siteDict[site]['file'] = [fileName]
                    siteDict[site]['seq'] = seq

siteFile = os.path.join(args.output, args.prefix + '.integrate.bed')
idStoreFile = os.path.join(args.output, args.prefix + '.idStore.txt')
loadFile = os.path.join(args.output, args.prefix + '.loadFile.txt')

with open(loadFile, 'w') as out:
    for key in keyList:
        for file in sorted(inputDict[key]):
            fileName = os.path.basename(file)
            row = [key, fileName]
            out.write('\t'.join(row) + '\n')

siteList = sorted(siteDict.keys())
with open(siteFile, 'w') as sf, open(idStoreFile, 'w') as df:
    sfRow = ["#Chro", "Start", "End", "SiteId", "MaxReadCount", "Strand", 
        "Seq", "SiteType", "RepNum", "Replicates", "PipelineNum", "Pipeline"]
    sf.write('\t'.join(sfRow) + '\n')
    for i in range(len(siteList)):
        site = siteList[i]
        siteSeq = siteDict[site]['seq']
        siteRepList = sorted(list(set(siteDict[site]['rep'])))
        siteTypeList = list(set(siteDict[site]['type']))
        siteRepNum = str(len(siteRepList))
        if len(siteTypeList) > 1:
            siteType = 'CIMS|CITS'
        else:
            siteType = siteTypeList[0]
        siteId = '{0}_{1}_{2}'.format(siteSeq,siteType.replace('|', '_'), i+1)
        siteRep = ','.join(siteRepList)
        keyList = sorted(list(set(siteDict[site]['key'])))
        key = ','.join(keyList)
        keyNum = str(len(keyList))
        readCount = str(siteDict[site]['read'])
        idCol = '|'.join(siteDict[site]['store'])
        fileName = '|'.join(siteDict[site]['file'])
        chro, start, end, strand = site.split('\t')
        sfRow = [chro, start, end, siteId, readCount, strand, siteSeq, siteType, siteRepNum, siteRep, keyNum, key]
        dfRow = [siteId, idCol, fileName]
        sf.write('\t'.join(sfRow) + '\n')
        df.write('\t'.join(dfRow) + '\n')
