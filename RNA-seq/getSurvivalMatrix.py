#!/usr/bin/env python3
import os
import sys
import argparse
import re
import numpy as np
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input', action='store', type=str, required=True,
                    help='input gene expression matrix from xena')
parser.add_argument('--clinc', action='store', type=str, required=True,
                    help='input clinical matrix')
parser.add_argument('--name', nargs='+', type=str,
                    help='gene names')
parser.add_argument('--id', nargs='+', type=str,
                    help='gene ids')
parser.add_argument('--upper', action='store', type=int,
                    default=75, help='percentile for high expression')
parser.add_argument('--lower', action='store', type=int,
                    default=25, help='percentile for low expression')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='output result')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

if bool(args.id) is False and bool(args.name) is False :
        parser.print_help()
        parser.exit()

#read expssion data
expData = pd.read_csv(args.input, sep='\t', header=0, index_col=0)

# filter with gene id
if bool(args.id):
    selectData = expData.loc[args.id]
else:
    selectData = expData[expData.GeneName.isin(args.name)]
    selectData = selectData.reset_index(drop=True)
    selectData = selectData.set_index("GeneName")

# delete index name and stack()
selectData.index.name = None
selectData = selectData.stack().unstack(0)

## get data with tumor samples
selectData = selectData.filter(regex='-\w+-0[1-9]\w', axis=0)

# read clinical data
clincData = pd.read_csv(args.clinc, sep='\t', header=0, index_col=0)
selectData.index.name = None
## get data with tumor samples
clincData = clincData.filter(regex='-\w+-0[1-9]\w', axis=0)

# common index between 2 dataframe
commonIndex = selectData.index.intersection(clincData.index)

selectData = selectData[selectData.index.isin(commonIndex)]
clincData = clincData[clincData.index.isin(commonIndex)]

# merge expression and clinical data 
mergeData = pd.merge(selectData, clincData, left_index=True, right_index=True, how='outer')

highExpList = list()
lowExpList = list()

if bool(args.id):
    geneList = args.id
else:
    geneList = args.gene

for gene in geneList:
    highExpData = mergeData[mergeData[gene] > np.percentile(mergeData[gene], args.upper)]
    lowExpData = mergeData[mergeData[gene] < np.percentile(mergeData[gene], args.lower)]
    highExpList.append(highExpData)
    lowExpList.append(lowExpData)

commonHighIndex = highExpList[0].index
commonLowIndex = lowExpList[0].index
high = highExpList[0]
low = lowExpList[0]
for i in range(1, len(highExpList)):
    commonHighIndex = high.index.intersection(highExpList[i].index)
    high = high[high.index.isin(commonHighIndex)]
    commonLowIndex = low.index.intersection(lowExpList[i].index)
    low = low[low.index.isin(commonLowIndex)]

high = high.assign(Category = ['high' for i in range(len(high))])
low = low.assign(Category = ['low' for i in range(len(low))])
contExpData = pd.concat([high, low], axis=0)

csv = contExpData.to_csv(index=True, index_label="Sample", sep="\t")
with open(args.output, 'w') as out:
    out.write(csv)
