#!/usr/bin/env python3
import os
import sys
import re
import numpy as np
from scipy import stats
import pandas as pd
import argparse
from multiprocessing import Pool, Manager
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--cole', action='store', type=int, required=True,
                    help='the end index of column for expression data (-1 indicate the end of the dataframe)')
parser.add_argument('--cols', action='store', type=int, required=True,
                    help='the start index of column for expression data')
parser.add_argument('--gene', action='store', nargs='+', type=str, required=True,
                    help='based gene list used for testing the co-expression network (if set as "all", then run program for all genes)')
parser.add_argument('--index', action='store', type=int, required=True,
                    help='the index of column used for pandas index (geneName, geneId, etc.)')
parser.add_argument('--contain', action='store', type=str,
                    help='regex for filtering the columns (included)')
parser.add_argument('--filter', action='store', type=str,
                    help='regex for filtering the columns (not included)')
parser.add_argument('--input', action='store', type=str, required=True,
                    help='input gene expression matrix (column:sample, row:gene)')
parser.add_argument('--mins', action='store', type=int,
                    default=10,
                    help='input gene expression matrix (column:sample, row:gene)')
parser.add_argument('--sep', action='store', type=str,
                    default='\t',
                    help='delimiter of columns to use')
parser.add_argument('--threads', action='store', type=int,
                    default=1,
                    help='threads to run the program')
parser.add_argument('--log2', action='store_true',
                    default=False,
                    help='used log2 to transform data')
parser.add_argument('--transpose', action='store_true',
                    default=False,
                    help='transpose the input matrix data')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='output result matrix')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def CallCoExpNetwork(data, igene, tgeneList, mins, cols, cole):
    # get expression data of input gene
    if cole == -1:
        igData = data.iloc[data.index == igene].iloc[:,cols:]
    else:
        igData = data.iloc[data.index == igene].iloc[:,cols:cole]
    # filter NA values
    igData = igData[igData.columns[~igData.isnull().all()]]
    
    coefList = list()
    for tgene in tgeneList:
        if tgene != igene:
            # get expression data of testing gene
            if cole == -1:
                geneData = data.iloc[data.index == tgene].iloc[:,cols:]
            else:
                geneData = data.iloc[data.index == tgene].iloc[:,cols:cole]
            coef = 0
            pval = 1
            sampleNum = 0
            for i in range(len(geneData)):
                # geneData may contain multiple rows
                tgData = geneData.iloc[[i]]
                # filter NA values
                tgData = tgData[tgData.columns[~tgData.isnull().all()]]
                ## get data from common columns and flatten 
                nparrA = igData[igData.columns & tgData.columns].to_numpy()[0]
                nparrB = tgData[igData.columns & tgData.columns].to_numpy()[0]
                tsampleNum = len(nparrA)
                if tsampleNum >= mins:
                    if args.log2 is True:
                        nparrA = np.log2(nparrA + 0.01)
                        nparrB = np.log2(nparrB + 0.01)
                    tcoef, tpvalue = stats.pearsonr(nparrA, nparrB)
                    if abs(coef) < abs(tcoef):
                        coef = tcoef
                        pvalue = tpvalue
                        sampleNum = tsampleNum
            if coef != 0:
                coefRow = [igene, tgene, str(coef), str(pvalue), str(sampleNum)]
                coefList.append(coefRow)
    return coefList

data = pd.read_csv(args.input, sep=args.sep, header=0, index_col=args.index)

# transpose data matrix if needed
if args.transpose is True:
    data = data.T

# filter column if needed
colNames = list(data.columns)

if bool(args.contain):
    colNames = list(filter(lambda x:bool(re.search(r'{0}'.format(args.contain), x)) is True, colNames ))
    data = data.filter(items=colNames)

if bool(args.filter):
    colNames = list(filter(lambda x:bool(re.search(r'{0}'.format(args.filter), x)) is False, colNames ))
    data = data.filter(items=colNames)

indexList = sorted(set(data.index.values), key=lambda x:str(x))
# run coexpression network
pool = Pool(processes=args.threads)
resultList = []
if args.gene == 'all':
    for i in range(len(indexList) - 1):
        igene = indexList[i]
        tgeneList = indexList[i+1:]
        result = pool.apply_async(CallCoExpNetwork, args=(data, igene, tgeneList, args.mins, args.cols, args.cole, ))
        resultList.append(result)
else:
    for igene in args.gene:
        if igene in indexList:
            tgeneList = list(filter(lambda x:x != args.gene, indexList))
            result = pool.apply_async(CallCoExpNetwork, args=(data, igene, tgeneList, args.mins, args.cols, args.cole, ))
            resultList.append(result)
        else:
            sys.error.write('No such gene ({0}) found in the expression matrix!'.format(igene))
            sys.exit()
pool.close()
pool.join()

with open(args.output, 'w') as out:
    row = ['inputGene', 'testGene', 'PCC', 'pvalue', 'sampleNum']
    out.write('\t'.join(row) + '\n')
    for result in resultList:
        coefList = result.get()
        for coefRow in coefList:
            out.write('\t'.join(coefRow) + '\n')
