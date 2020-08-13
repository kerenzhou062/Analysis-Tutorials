#!/usr/bin/env python3
import os
import sys
import re
import numpy as np
from scipy import stats
import pandas as pd
import argparse
from multiprocessing import Pool
import time

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--cole', action='store', type=int,
                    default=0,
                    help='the end index of column of DataFrame (except for the index column, -1 indicate the end of the dataframe)')
parser.add_argument('--cols', action='store', type=int,
                    default=-1,
                    help='the start index of column of DataFrame (except for the index column)')
parser.add_argument('--gene', action='store', nargs='+', type=str, required=True,
                    help='based gene list used for testing the co-expression network (if set as "all", then run program for all genes)')
parser.add_argument('--index', action='store', type=int, required=True,
                    help='the index of column of input matrix used for DataFrame index (geneName, geneId, etc.)')
parser.add_argument('--contain', action='store', type=str,
                    help='regex for filtering the columns (included)')
parser.add_argument('--filter', action='store', type=str,
                    help='regex for filtering the columns (excluded)')
parser.add_argument('--input', action='store', type=str, required=True,
                    help='input gene expression matrix (column:sample, row:gene)')
parser.add_argument('--minsize', action='store', type=int,
                    default=10,
                    help='input gene expression matrix (column:sample, row:gene)')
parser.add_argument('--operator', action='store', nargs='*', type=str,
                    choices=['>', '>=', '<', '<=', '!='],
                    help='operators used for filtering the DataFrame')
parser.add_argument('--opertype', action='store', type=str,
                    default='multiple',
                    choices=['single', 'multiple'],
                    help='apply operators on single row or multiple row of DataFrame')
parser.add_argument('--operval', action='store', nargs='*', type=float,
                    help='filtering values for operator (corresponding to --operator)')
parser.add_argument('--sep', action='store', type=str,
                    default='\t',
                    help='delimiter of columns to use')
parser.add_argument('--threads', action='store', type=int,
                    default=1,
                    help='threads to run the program')
parser.add_argument('--log2', action='store_true',
                    default=False,
                    help='used log2 to transform data')
parser.add_argument('--time', action='store_true',
                    default=False,
                    help='report running time')
parser.add_argument('--transpose', action='store_true',
                    default=False,
                    help='transpose the input matrix data before running program')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='output result matrix')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def FilterMatrix(df, operator, operval):
    # remove columns with na
    df = df.dropna(axis=1, how='any')
    ## filter columns by row values
    if bool(operator) and bool(operval) :
        if operator == '<':
            df = df[df.columns[(df < operval).any()]]
        elif operator == '<=':
            df = df[df.columns[(df <= operval).any()]]
        elif operator == '>':
            df = df[df.columns[(df > operval).any()]]
        elif operator == '>=':
            df = df[df.columns[(df >= operval).any()]]
        elif operator == '!=':
            df = df[df.columns[(df != operval).any()]]
    return df

def CallCoExpNetwork(data, igene, tgeneList, minSize, operators, opervals, opertype):
    # get expression data of input gene
    igAllData = data.iloc[data.index == igene]
    # filter out values
    if opertype == 'single':
        for i in range(len(operators)):
            igAllData = FilterMatrix(igAllData, operators[i], opervals[i])

    coefList = list()
    # igAllData may contain multiple rows
    for i in range(len(igAllData)):
        igData = igAllData.iloc[[i]]
        for tgene in tgeneList:
            if tgene != igene:
                # get expression data of testing gene
                tgAllData = data.iloc[data.index == tgene]
                rho = 0
                pval = 1
                sampleSize = 0
                for j in range(len(tgAllData)):
                    # tgAllData may contain multiple rows
                    tgData = tgAllData.iloc[[j]]
                    # filter out values
                    if opertype == 'single':
                        for k in range(len(operators)):
                            tgData = FilterMatrix(tgData, operators[k], opervals[k])
                            ## get data from common columns and flatten 
                            nparrA = igData[igData.columns & tgData.columns].to_numpy()[0]
                            nparrB = tgData[igData.columns & tgData.columns].to_numpy()[0]
                    else:
                        # Concatenate igData and tgData
                        mergeData = pd.concat([igData, tgData], axis=0)
                        for k in range(len(operators)):
                            # filter out values on multiple rows
                            mergeData = FilterMatrix(mergeData, operators[k], opervals[k])
                            ## get data from common columns and flatten 
                            nparrA = mergeData.iloc[[0]].to_numpy()[0]
                            nparrB = mergeData.iloc[[1]].to_numpy()[0]
                    ## to avoid 0 elements
                    if len(nparrA) == 0 or len(nparrB) == 0:
                        continue
                    ## to avoid PearsonRConstantInputWarning: constant values
                    if np.all(nparrA == nparrA[0]) or np.all(nparrB == nparrB[0]):
                        continue
                    tsampleSize = len(nparrA)
                    if tsampleSize >= minSize:
                        if args.log2 is True:
                            nparrA = np.log2(nparrA + 0.01)
                            nparrB = np.log2(nparrB + 0.01)
                        trho, tpvalue = stats.pearsonr(nparrA, nparrB)
                        if abs(rho) < abs(trho):
                            rho = trho
                            pvalue = tpvalue
                            sampleSize = tsampleSize
                if rho != 0:
                    coefRow = [igene, tgene, str(rho), str(pvalue), str(sampleSize)]
                    coefList.append(coefRow)
    return coefList

# record time
stime = time.time()
# check --operator and --operval

if bool(args.operator) and bool(args.operval):
    if len(args.operator) != len(args.operval):
        sys.stderr.write('Errors in --operator and --operval!\n')
        sys.exit()
else:
    args.operator = []
    args.operval = []

# read input data into a matrix
data = pd.read_csv(args.input, sep=args.sep, header=0, index_col=args.index)

# transpose data matrix if needed
if args.transpose is True:
    data = data.T

# get expression data from column-cols to column-cole
if args.cole == -1:
    data = data.iloc[:,args.cols:]
else:
    data = data.iloc[:,args.cols:args.cole]

# filter column if needed
colNames = list(data.columns)

if bool(args.contain):
    colNames = list(filter(lambda x:bool(re.search(r'{0}'.format(args.contain), x)) is True, colNames ))
    data = data.filter(items=colNames)

if bool(args.filter):
    colNames = list(filter(lambda x:bool(re.search(r'{0}'.format(args.filter), x)) is False, colNames ))
    data = data.filter(items=colNames)

# get gene list
indexList = sorted(set(data.index.values), key=lambda x:str(x))

# run coexpression network
pool = Pool(processes=args.threads)
resultList = []
if args.gene == 'all':
    for i in range(len(indexList) - 1):
        igene = indexList[i]
        tgeneList = indexList[i+1:]
        result = pool.apply_async(CallCoExpNetwork, args=(data, igene, tgeneList, args.minsize, args.operator, args.operval, args.opertype,))
        resultList.append(result)
else:
    for igene in args.gene:
        if igene not in indexList:
            sys.error.write('No such gene ({0}) found in the expression matrix!\n'.format(igene))
            sys.exit()
    for igene in args.gene:
        tgeneList = list(filter(lambda x:x != args.gene, indexList))
        result = pool.apply_async(CallCoExpNetwork, args=(data, igene, tgeneList, args.minsize, args.operator, args.operval, args.opertype,))
        resultList.append(result)
pool.close()
pool.join()

with open(args.output, 'w') as out:
    row = ['inputGene', 'testGene', 'pcc', 'pvalue', 'sampleSize']
    out.write('\t'.join(row) + '\n')
    for result in resultList:
        coefList = result.get()
        for coefRow in coefList:
            out.write('\t'.join(coefRow) + '\n')

if args.time:
    print("Running time is: %s seconds." % (time.time() - stime))
