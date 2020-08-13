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
parser.add_argument('--cole', action='store', type=int,
                    default=0,
                    help='the end index of column for expression data (except for the index column, -1 indicate the end of the dataframe)')
parser.add_argument('--cols', action='store', type=int,
                    default=-1,
                    help='the start index of column for expression matrix (except for the index column)')
parser.add_argument('--gene', action='store', nargs='+', type=str, required=True,
                    help='based gene list used for testing the co-expression network (if set as "all", then run program for all genes)')
parser.add_argument('--index', action='store', type=int, required=True,
                    help='the index of column used for pandas index (geneName, geneId, etc.)')
parser.add_argument('--contain', action='store', type=str,
                    help='regex for filtering the columns (included)')
parser.add_argument('--filter', action='store', type=str,
                    help='regex for filtering the columns (excluded)')
parser.add_argument('--input', action='store', type=str, required=True,
                    help='input gene expression matrix (column:sample, row:gene)')
parser.add_argument('--mins', action='store', type=int,
                    default=10,
                    help='input gene expression matrix (column:sample, row:gene)')
parser.add_argument('--operator1', action='store', type=str,
                    choices=['>', '>=', '<', '<=', '!='],
                    help='operator used for filtering the matrix data')
parser.add_argument('--operval1', action='store', type=float,
                    help='filtering value for operator')
parser.add_argument('--operator2', action='store', type=str,
                    choices=['>', '>=', '<', '<=', '!='],
                    help='operator used for filtering the matrix data')
parser.add_argument('--operval2', action='store', type=float,
                    help='filtering value for operator')
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

def FilterMatrix(df, operator, operval):
    df = df[df.columns[~df.isnull().all()]]
    ## filter columns by row values
    if bool(operator) and bool(operval) :
        if operator == '<':
            df = df[df.columns[df.iloc[0, :] < operval]]
        elif operator == '<=':
            df = df[df.columns[df.iloc[0, :] <= operval]]
        elif operator == '>':
            df = df[df.columns[df.iloc[0, :] > operval]]
        elif operator == '>=':
            df = df[df.columns[df.iloc[0, :] >= operval]]
        elif operator == '!=':
            df = df[df.columns[df.iloc[0, :] != operval]]
    return df

def CallCoExpNetwork(data, igene, tgeneList, mins, operators, opervals):
    # get expression data of input gene
    igData = data.iloc[data.index == igene]
    # filter NA values
    igData = FilterMatrix(igData, operators[0], opervals[0])
    igData = FilterMatrix(igData, operators[1], opervals[1])

    coefList = list()
    for tgene in tgeneList:
        if tgene != igene:
            # get expression data of testing gene
            geneData = data.iloc[data.index == tgene]
            coef = 0
            pval = 1
            sampleNum = 0
            for i in range(len(geneData)):
                # geneData may contain multiple rows
                tgData = geneData.iloc[[i]]
                # filter NA values
                tgData = FilterMatrix(tgData, operators[0], opervals[0])
                tgData = FilterMatrix(tgData, operators[1], opervals[1])
                ## get data from common columns and flatten 
                nparrA = igData[igData.columns & tgData.columns].to_numpy()[0]
                nparrB = tgData[igData.columns & tgData.columns].to_numpy()[0]
                ## to avoid 0 elements
                if len(nparrA) == 0 or len(nparrB) == 0:
                    continue
                ## to avoid PearsonRConstantInputWarning: constant values
                if np.all(nparrA == nparrA[0]) or np.all(nparrB == nparrB[0]):
                    continue
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

## read input data into a matrix
data = pd.read_csv(args.input, sep=args.sep, header=0, index_col=args.index)

# transpose data matrix if needed
if args.transpose is True:
    data = data.T

## get expression data from column-cols to column-cole
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

indexList = sorted(set(data.index.values), key=lambda x:str(x))
# run coexpression network
operators = [args.operator1, args.operator2]
opervals = [args.operval1, args.operval2]

pool = Pool(processes=args.threads)
resultList = []
if args.gene == 'all':
    for i in range(len(indexList) - 1):
        igene = indexList[i]
        tgeneList = indexList[i+1:]
        result = pool.apply_async(CallCoExpNetwork, args=(data, igene, tgeneList, args.mins, operators, opervals, ))
        resultList.append(result)
else:
    for igene in args.gene:
        if igene in indexList:
            tgeneList = list(filter(lambda x:x != args.gene, indexList))
            result = pool.apply_async(CallCoExpNetwork, args=(data, igene, tgeneList, args.mins, operators, opervals, ))
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
