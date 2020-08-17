#!/usr/bin/env python3
import os
import sys
import re
import argparse
import numpy as np
from scipy import stats
import pandas as pd
from multiprocessing import Pool
from collections import defaultdict
import time

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--faxis', action='store', type=int,
                    default=0,
                    choices=[0, 1, 2],
                    help='axis to apply the --filter (0:row, 1:column, 2:both)')
parser.add_argument('--fillna', action='store', type=float,
                    help='fill the "nan" values with # in the DataFrame')
parser.add_argument('--filter', action='store', type=str,
                    help='keep labels from axis for which re.search(regex, label) == True (After transpose if needed)')
parser.add_argument('--gene', action='store', nargs='+', type=str, required=True,
                    help='based gene list used for testing the co-expression network (if set as "all", then run program for all genes)')
parser.add_argument('--index', action='store', type=int, required=True,
                    help='used #th column (0-based) of --input for DataFrame index (geneName, geneId, etc.)')
parser.add_argument('--input', action='store', type=str, required=True,
                    help='input gene expression matrix (column:sample, row:gene, "na" values will be filtered before calculation)')
parser.add_argument('--method', action='store', type=str,
                    choices=['spearman', 'pearson'],
                    default='pearson',
                    help='method to calculate the coefficient')
parser.add_argument('--minsize', action='store', type=int,
                    default=10,
                    help='minimun sample size for calculating pearson correlation coefficient')
parser.add_argument('--operator', action='store', nargs='*', type=str,
                    choices=['>', '>=', '<', '<=', '!='],
                    help='operators used for filtering the DataFrame')
parser.add_argument('--operval', action='store', nargs='*', type=float,
                    help='filtering values for operator (corresponding to --operator)')
parser.add_argument('--opertype', action='store', type=str,
                    default='paired',
                    choices=['single', 'paired'],
                    help='apply operators on single row or paired row of DataFrame')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='output result matrix')
parser.add_argument('--sep', action='store', type=str,
                    default='\t',
                    help='delimiter of columns to use')
parser.add_argument('--slice', action='store', nargs=2, type=int,
                    default=[0, -1],
                    help='slice the DataFrame from 1#th to 2#th row (-1 indicate the end of the dataframe)')
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
                    help='transpose the DataFrame to fit --input data structure before calculating pcc')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def FiltDataframe(df, operator, operval, opertype):
    ## filter columns by row values
    if (type(operator) is str) and (type(operval) is float):
        if opertype == 'single':
            if operator == '<':
                df =df[df.iloc[:,[0]] < operval]
            elif operator == '<=':
                df = df[df.iloc[:,[0]] <= operval]
            elif operator == '>':
                df =df[df.iloc[:,[0]] > operval]
            elif operator == '>=':
                df = df[df.iloc[:,[0]] >= operval]
            elif operator == '!=':
                df = df[df.iloc[:,[0]] != operval]
        else:
            cols = df.columns
            if operator == '<':
                df = df[(df[cols[0]] < operval) | (df[cols[1]] < operval)]
            elif operator == '<=':
                df = df[(df[cols[0]] <= operval) | (df[cols[1]] <= operval)]
            elif operator == '>':
                df = df[(df[cols[0]] > operval) | (df[cols[1]] > operval)]
            elif operator == '>=':
                df = df[(df[cols[0]] >= operval) | (df[cols[1]] >= operval)]
            elif operator == '!=':
                df = df[(df[cols[0]] != operval) | (df[cols[1]] != operval)]
    # remove columns with na
    df = df.dropna(axis=0, how='any')
    return df

def CallCoExpNetwork(data, qgene, tgeneList, minSize, operators, opervals, opertype, method):
    # get expression data of input gene
    igAllData = data[[qgene]]
    coefDict = defaultdict(list)
    # igAllData may contain multiple rows
    for i in range(len(igAllData.columns)):
        igData = igAllData.iloc[:, [i]]
        # filter out values
        if opertype == 'single':
            for i in range(len(operators)):
                igData = FiltDataframe(igData, operators[i], opervals[i], opertype)
        for tgene in tgeneList:
            rho = 0
            pval = 1
            sampleSize = 0
            # get expression data of testing gene
            tgAllData = data[[tgene]]
            for j in range(len(tgAllData.columns)):
                # tgAllData may contain multiple rows
                tgData = tgAllData.iloc[:, [j]]
                # filter out values
                if opertype == 'single':
                    for k in range(len(operators)):
                        tgData = FiltDataframe(tgData, operators[k], opervals[k], opertype)
                        ## get data from common columns and flatten
                    mergeData = pd.concat([igData, tgData], axis=1)
                    mergeData = FiltDataframe(mergeData, False, False, False)
                else:
                    # Concatenate igData and tgData
                    mergeData = pd.concat([igData, tgData], axis=1)
                    mergeData = FiltDataframe(mergeData, False, False, False)
                    for k in range(len(operators)):
                        # filter out values on multiple rows
                        mergeData = FiltDataframe(mergeData, operators[k], opervals[k], opertype)
                nparrA = mergeData[qgene].astype('float64').to_numpy()
                nparrB = mergeData[tgene].astype('float64').to_numpy()
                ## to avoid 0 elements
                if len(nparrA) == 0:
                    continue
                ## to avoid PearsonRConstantInputWarning: constant values
                if np.all(nparrA == nparrA[0]) or np.all(nparrB == nparrB[0]):
                    continue
                tsampleSize = len(nparrA)
                if tsampleSize >= minSize:
                    if args.log2 is True:
                        nparrA = np.log2(nparrA + 0.01)
                        nparrB = np.log2(nparrB + 0.01)
                    if method == 'pearson':
                        trho, tpvalue = stats.pearsonr(nparrA, nparrB)
                    elif method == 'spearman':
                        trho, tpvalue = stats.spearmanr(nparrA, nparrB)
                    if abs(rho) < abs(trho):
                        rho = trho
                        pvalue = tpvalue
                        sampleSize = tsampleSize
                if rho != 0:
                    if tgene not in coefDict:
                        coefDict[tgene] = [rho, pvalue, sampleSize]
                    else:
                        if abs(rho) > coefDict[tgene][2]:
                            coefDict[tgene] = [rho, pvalue, sampleSize]
    return [qgene, coefDict]

# record start time
stime = time.time()

# check --operator and --operval
if args.operator is not None and args.operval is not None:
    if len(args.operator) != len(args.operval):
        sys.stderr.write('Errors in --operator and --operval!\n')
        sys.exit()
else:
    args.operator = []
    args.operval = []

# read input data into a matrix
data = pd.read_csv(args.input, sep=args.sep, header=0, index_col=args.index, engine='c', memory_map=True, converters={args.index:str})

# drop missing values in index
data = data.loc[data.index.dropna()]

# transpose data matrix if needed
if args.transpose is True:
    data = data.T
    # drop missing values in index after transpose
    data = data.loc[data.index.dropna()]

# fill nan value with --fillna
if args.fillna is not None:
    data = data.fillna(args.fillna)

# get expression data from row-start to row-end
if args.slice[0] < 0:
    args.slice[0] = 0

if args.slice[1] == -1 or args.slice[1] > len(data):
    args.slice[1] = len(data) - 1

data = data.iloc[args.slice[0]:args.slice[1], :]

# filter column and index if needed
if bool(args.filter):
    if args.faxis == 2:
        data = data.filter(regex=args.filter, axis=0)
        data = data.filter(regex=args.filter, axis=1)
    else:
        data = data.filter(regex=args.filter, axis=args.faxis)

# get final gene list and remove the duplicate records
geneList = sorted(set(data.columns))

# run coexpression network
pool = Pool(processes=args.threads)
resultList = []

if 'all' in args.gene:
    for i in range(len(geneList) - 1):
        qgene = geneList[i]
        tgeneList = geneList[i+1:]
        result = pool.apply_async(CallCoExpNetwork, args=(data, qgene, tgeneList, args.minsize, args.operator, args.operval, args.opertype, args.method))
        resultList.append(result)
else:
    for qgene in args.gene:
        if qgene not in geneList:
            sys.stderr.write('No such gene ({0}) found in the expression matrix!\n'.format(qgene))
            sys.exit()
    for qgene in args.gene:
        tgeneList = list(filter(lambda x:x not in args.gene, geneList))
        result = pool.apply_async(CallCoExpNetwork, args=(data, qgene, tgeneList, args.minsize, args.operator, args.operval, args.opertype, args.method))
        resultList.append(result)
pool.close()
pool.join()

with open(args.output, 'w') as out:
    row = ['inputGene', 'testGene', 'pcc', 'pvalue', 'sampleSize']
    out.write('\t'.join(row) + '\n')
    for result in resultList:
        qgene, coefDict = result.get()
        for tgene in sorted(coefDict.keys()):
            coefRow = [qgene, tgene] + coefDict[tgene]
            out.write('\t'.join(map(str, coefRow)) + '\n')

if args.time:
    print("Running time is: %s seconds." % (time.time() - stime))
