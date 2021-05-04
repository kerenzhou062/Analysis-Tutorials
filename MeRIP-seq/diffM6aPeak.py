#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
import tempfile
import subprocess
import math
import numpy as np
import scipy.stats as stats
from multiprocessing import Pool
## used R packages
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-c', '--control', action='store', type=str,
                    required=True,
                    help='M6A peak.xls of control group from exomePeak')
parser.add_argument('-f', '--filter', action='store', type=int,
                    default=-1,
                    help='filter out genes less than # counts across 1/2 samples (--filter in DESeq2Gene.R)')
parser.add_argument('-l', '--library', action='store', type=str,
                    choices=['unstranded', 'fr-firstrand','fr-secondstrand'],
                    default='fr-firstrand',
                    help='Library protocols of bam')
parser.add_argument('-o', '--output', action='store', type=str,
                    required=True,
                    help='The output peak.custom.diff.xls')
parser.add_argument('-p', '--pool', action='store', type=str,
                    required=True,
                    help='M6A peak.xls of ([IP.treat + IP.control] vs [input.treat + input.control]) group from exomePeak')
parser.add_argument('-s', '--shrink', action='store', type=str, choices=['none', 'apeglm', 'ashr'],
                    default='none',
                    help='The method used for shrinking the fold change (lfcShrink() in DESeq2)')
parser.add_argument('-t', '--treat', action='store', type=str,
                    required=True,
                    help='M6A peak.xls of treatment group from exomePeak')
parser.add_argument('--bamcip', nargs='+', type=str,
                    required=True,
                    help='bam files of conrtol IP samples')
parser.add_argument('--bamcinput', nargs='+', type=str,
                    required=True,
                    help='bam files of conrtol input samples')
parser.add_argument('--bamtip', nargs='+', type=str,
                    required=True,
                    help='bam files of treated IP samples')
parser.add_argument('--bamtinput', nargs='+', type=str,
                    required=True,
                    help='bam files of treated input samples')
parser.add_argument('--pairend', action='store_true', 
                    default=False,
                    help='the input bams are pair-end')
parser.add_argument('--expMtx', action='store', type=str,
                    help='The expression matrix from buildExpMatrix.py')
parser.add_argument('--degMtx', action='store', type=str,
                    help='The differentially expressed gene matrix from DESeq2Gene.R')
parser.add_argument('--cntKey', action='store', type=str,
                    help='Keyword for names of control samples (required by --expMtx)')
parser.add_argument('--trtKey', action='store', type=str,
                    help='Keyword for names of treatment samples (required by --expMtx)')
parser.add_argument('--thread', action='store', type=int,
                    default=1,
                    help='The number of threads used to run in parallel (will greatly reduce running time)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def GetPeakReadInBam(mapArgsList):
    peakFile, bamFile, uniqBamIndex, strand, pairendFlag = mapArgsList
    ## get chromosome names from bam.bai
    gsizeTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
    command = 'samtools idxstats {}'.format(bamFile)
    idxstatsList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
    idxstatsList = list(filter(lambda x: (bool(x) and x != '*'), idxstatsList))
    totalMapReadNum = 0
    with open(gsizeTmp.name, 'w') as temp:
        for line in idxstatsList:
            if bool(line) is False:
                continue
            row = line.strip().split('\t')
            chrom = row[0]
            if re.match(r'ERCC', chrom):
                continue
            chromLen = row[1]
            totalMapReadNum += int(row[2])
            temp.write('\t'.join([chrom, chromLen]) + '\n')
    ## get final peak-read dictionary
    prDict = {}
    if pairendFlag is True:
        ## sort peaks by chromosomes order in bam
        sortedPeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
        command = 'bedtools sort -i {} -g {} > {}'.format(peakFile, gsizeTmp.name, sortedPeakTmp.name)
        __ = subprocess.check_output(command, shell=True)
        ## intersect peaks with sorted bam
        if strand == '-s':
            reverseFlag = False
            unstrandFlag = False
        elif strand == '-S':
            reverseFlag = True
            unstrandFlag = False
        else:
            reverseFlag = False
            unstrandFlag = True
        command = 'bedtools intersect -abam {} -b {} -g {} -sorted -split -bed -wa -wb'.format(bamFile, sortedPeakTmp.name, gsizeTmp.name)
        intersecctList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
        sortedPeakTmp.close()
        gsizeTmp.close()
        ## construct peakid-bam-reads
        intersectDict = defaultdict(dict)
        for line in intersecctList:
            if bool(line) is False:
                continue
            row = line.strip().split('\t')
            readid = row[3]
            readStrand = row[5]
            peakid = row[15]
            peakStrand = row[17]
            readName = readid[:-2]
            mate = readid[-2:]
            if reverseFlag is True:
                if (mate == '/1' and readStrand == peakStrand) or (mate == '/2' and readStrand != peakStrand):
                    continue
            elif reverseFlag is False and unstrandFlag is False:
                if (mate == '/1' and readStrand != peakStrand) or (mate == '/2' and readStrand == peakStrand):
                    continue
            readid = readName
            intersectDict[peakid][readName] = 1
        for peakid in sorted(intersectDict.keys()):
            prDict[peakid] = len(intersectDict[peakid].keys())
    else:
        ## sort peaks by chromosome && intersect peaks with sorted bam
        command = 'bedtools sort -i {} -g {} | '.format(peakFile, gsizeTmp.name)
        command += 'bedtools intersect -a stdin -b {} -g {} -sorted -split -c {}'.format(bamFile, gsizeTmp.name, strand)
        peakReadsList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
        gsizeTmp.close()
        for line in peakReadsList:
            if bool(line) is False:
                continue
            row = line.strip().split('\t')
            peakid = row[3]
            readsNum = row[-1]
            prDict[peakid] = int(readsNum)
    return [uniqBamIndex, prDict]

def GetAjustPvalBH(pvalList):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    orderDescend = p.argsort()[::-1]
    orderOrig = orderDescend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[orderDescend]))
    return q[orderOrig]

def RunFisherExactTest(peakid):
    controlList = [np.mean(norCountsDict[peakid]['control']['IP']), np.mean(norCountsDict[peakid]['control']['input'])]
    treatList = [np.mean(norCountsDict[peakid]['treat']['IP']), np.mean(norCountsDict[peakid]['treat']['input'])]
    oddsratio, pvalue = stats.fisher_exact([controlList, treatList])
    return [peakid, pvalue]

args.output = os.path.realpath(args.output)
thread = args.thread
pairendFlag = args.pairend
rstats = importr('stats')

##build bam dict
ipTypeList  = ['IP', 'input', 'IP', 'input']
treatTypeList = ['control', 'control', 'treat', 'treat']
bamList = [args.bamcip, args.bamcinput, args.bamtip, args.bamtinput]

uniqBamIndexList = []
bamDict = defaultdict(dict)
bamDict['bamToIndex'] = defaultdict(dict)
bamDict['indexToInfo'] = defaultdict(dict)
for i in range(4):
    ipType = ipTypeList[i]
    treatType = treatTypeList[i]
    bamFileList = bamList[i]
    for j in range(len(bamFileList)):
        bamFile = bamFileList[j]
        uniqBamIndex = '_'.join([str(i), str(j)])
        if treatType not in bamDict:
            bamDict[treatType] = defaultdict(list)
        bamDict['bamToIndex'][bamFile] = uniqBamIndex
        bamDict['indexToInfo'][uniqBamIndex]['ipType'] = ipType
        bamDict['indexToInfo'][uniqBamIndex]['treatType'] = treatType
        uniqBamIndexList.append(uniqBamIndex)

nameRow = list()
controlPeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
with open(args.control, 'r') as f, open(controlPeakTmp.name, 'w') as temp:
    nameRow = f.readline().strip().split('\t')
    nameRow[-3] = 'peak.pval'
    nameRow[-1] = 'peak.fold.enrichment.exomePeak'
    nameRow.append("peak.source")
    for line in f:
        temp.write(line)

treatPeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
with open(args.treat, 'r') as f, open(treatPeakTmp.name, 'w') as temp:
    __ = f.readline()
    for line in f:
        temp.write(line)

## get peaks only found in control or treated groups
command = 'bedtools intersect -a {} -b {} -f 0.8 -F 0.8 -e -v -s -split'.format(controlPeakTmp.name, treatPeakTmp.name)
controlUniqPeakList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')

command = 'bedtools intersect -a {} -b {} -f 0.8 -F 0.8 -e -v -s -split'.format(treatPeakTmp.name, controlPeakTmp.name)
treatUniqPeakList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')

controlPeakTmp.close()
treatPeakTmp.close()

# identify uniquely represented peaks in control or treat groups
sourceRow = []
combineRow = []
for line in controlUniqPeakList:
    if bool(line) is False:
        continue
    row = line.strip().split('\t')
    sourceRow.append('control')
    combineRow.append(row)

for line in treatUniqPeakList:
    if bool(line) is False:
        continue
    row = line.strip().split('\t')
    sourceRow.append('treated')
    combineRow.append(row)

# combine with pool peak results
combinePeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
poolPeakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
with open(combinePeakTmp.name, 'w') as temp:
    for row in combineRow:
        temp.write('\t'.join(row) + '\n')

with open(args.pool, 'r') as f, open(poolPeakTmp.name, 'w') as temp:
    __ = f.readline()
    for line in f:
        temp.write(line)

# remove peaks that overlap with peaks that only found in control or treated groups
command = 'bedtools intersect -a {} -b {} -f 0.8 -F 0.8 -e -v -s -split'.format(poolPeakTmp.name, combinePeakTmp.name)
poolUniqPeakList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')

combinePeakTmp.close()
poolPeakTmp.close()

for line in poolUniqPeakList:
    if bool(line) is False:
        continue
    row = line.strip().split('\t')
    sourceRow.append('pool')
    combineRow.append(row)

# build expression dict
expSampleDict = defaultdict(list)
expDict = defaultdict(dict)
if bool(args.expMtx):
    with open(args.expMtx, 'r') as f:
        expNameRow = list(filter(lambda x: bool(re.search(r'_CQV', x)) is False, f.readline().split('\t')))
        rowLength = len(expNameRow)
        for i in range(1,rowLength):
            if bool(args.cntKey) and bool(args.trtKey):
                if bool(re.search(r'{}'.format(args.cntKey), expNameRow[i])):
                    expSampleDict['control'].append(i)
                if bool(re.search(r'{}'.format(args.trtKey), expNameRow[i])):
                    expSampleDict['treat'].append(i)
            else:
                sys.stderr.write('--cntKey and --trtKey required by --expMtx!')
                sys.exit()
        for line in f:
            row = line.strip().split('\t')
            geneId = row[0].split('.')[0]
            expDict[geneId]['control'] = sum(map(lambda x: float(row[x]), expSampleDict['control'])) / len(expSampleDict['control'])
            expDict[geneId]['treat'] = sum(map(lambda x: float(row[x]), expSampleDict['treat'])) / len(expSampleDict['treat'])

# build differentlly expressed gene dict
geneDeDict = defaultdict(dict)
if bool(args.degMtx):
    with open(args.degMtx, 'r') as f:
        __ = f.readline()
        for line in f:
            row = line.strip().split('\t')
            geneId = row[0].split('.')[0]
            log2fc = '{0:.4f}'.format(float(row[-4]))
            pval = row[-2]
            if pval == 'NA':
                pval = '1'
            geneDeDict[geneId] = [log2fc, pval]

# prepare peak dictionary
peakidList = list()
peakDict = defaultdict(dict)
for i in range(len(combineRow)):
    geneId = combineRow[i][3].split('.')[0]
    combineRow[i][3] = '|'.join([combineRow[i][3], str(i)])
    peakid = combineRow[i][3]
    peakDict[peakid]['row'] = combineRow[i]
    peakDict[peakid]['source'] = sourceRow[i]
    peakidList.append(peakid)
    peakDict[peakid]['gene'] = geneId
    ## add expression
    if bool(expDict):
        expList = ['NA', 'NA']
        if 'control' in expDict[geneId]:
            expList[0] = '{0:.4f}'.format(expDict[geneId]['control'])
        if 'treat' in expDict[geneId]:
            expList[1] = '{0:.4f}'.format(expDict[geneId]['treat'])
        peakDict[peakid]['geneExp'] = expList
    else:
        peakDict[peakid]['geneDe'] = list()
    ## add differential expression
    if bool(geneDeDict):
        if geneId not in geneDeDict:
            peakDict[peakid]['geneDe'] = ['0', '1']
        else:
            peakDict[peakid]['geneDe'] = geneDeDict[geneId]
    else:
        peakDict[peakid]['geneDe'] = list()

outputDir = os.path.dirname(args.output)
genePeakKeptDict = defaultdict(dict)
if bool(args.bamcip):
    ## calculate reads from samples that cover peaks
    if args.library == 'unstranded':
        strand = ''
    elif args.library == 'fr-firstrand':
        strand = '-S'
    else:
        strand = '-s'
    peakTmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
    with open(peakTmp.name, 'w') as temp:
        for row in combineRow:
            temp.write('\t'.join(row) + '\n')
    # run peaks-reads in bams in parallel
    poolMapArgsList = []
    for i in range(len(bamList)):
        for bamFile in bamList[i]:
            uniqBamIndex = bamDict['bamToIndex'][bamFile]
            poolMapArgsList.append([peakTmp.name, bamFile, uniqBamIndex, strand, pairendFlag])
    peakReadsResultList = []
    with Pool(processes=thread) as pool:
        for i in pool.imap_unordered(GetPeakReadInBam, poolMapArgsList):
            peakReadsResultList.append(i)
    ## delete peak temp
    peakTmp.close()
    ## construct final peak-reads dictionary
    peakReadDict = defaultdict(dict)
    for result in peakReadsResultList:
        uniqBamIndex, prDict = result
        for peakid in sorted(prDict.keys()):
            peakReadDict[peakid][uniqBamIndex] = prDict[peakid]
    
    ## construct the sample matrix requried by DESeq2Gene.R
    sampleMtx = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
    with open(sampleMtx.name, 'w') as temp:
        row = ['sample', 'condition', 'antibody\n']
        temp.write('\t'.join(row))
        for uniqBamIndex in bamDict['indexToInfo'].keys():
            ipType = bamDict['indexToInfo'][uniqBamIndex]['ipType']
            treatType = bamDict['indexToInfo'][uniqBamIndex]['treatType']
            row = ["X" + uniqBamIndex, treatType, ipType + '\n']
            temp.write('\t'.join(row))
    ## construct the peak-counts matrix
    peakCountsMtx = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
    #peakCountsMtx = 'peak.counts.matrix'
    with open(peakCountsMtx.name, 'w') as temp:
        row = ['gene_id'] + list(map(lambda x: "X" + x, uniqBamIndexList))
        temp.write('\t'.join(row) + '\n')
        for peakid in sorted(peakReadDict.keys()):
            row = [peakid]
            for uniqBamIndex in uniqBamIndexList:
                if uniqBamIndex not in peakReadDict[peakid]:
                    counts = 0
                else:
                    counts = peakReadDict[peakid][uniqBamIndex]
                row.append(str(counts))
            temp.write('\t'.join(row) + '\n')
    ## run DESeq2Gene.R, m6a_peak_diff.DESeq2.txt, m6a_peak_diff.sig.DESeq2.txt, m6a_peak_diff.normalized.counts.txt
    command = 'DESeq2Gene.R --sampleMtx {} --counts {} \
        --formula "condition + antibody + condition:antibody" --relevel "condition:control,antibody:input" \
        --name "conditiontreat.antibodyIP" -f "{}" --norcounts \
        --test "Wald" --pval 1 --adjp 0.1 --shrink "{}" --prefix "m6a_peak_diff" --output {}'.format(sampleMtx.name, peakCountsMtx.name, args.filter, args.shrink, outputDir)
    __ = subprocess.check_output(command, shell=True)
    ## delete temporary files
    sampleMtx.close()
    peakCountsMtx.close()
    deseqResFile = os.path.join(outputDir, 'm6a_peak_diff.DESeq2.txt')
    with open(deseqResFile, 'r') as f:
        row = f.readline().strip().split('\t')
        headerDict = defaultdict(dict)
        for i in range(len(row)):
            headerDict[row[i]] = i
        for line in f:
            row = line.strip().split('\t')
            peakid = row[headerDict['uniqId']]
            baseMean = row[headerDict['baseMean']]
            if float(baseMean) == 0:
                continue
            log2Fold = row[headerDict['log2FoldChange']]
            #pval = row[headerDict['pvalue']]
            #pvalAdj = row[5]
            valueRow = [baseMean, log2Fold]
            peakDict[peakid]['peakDiffDeseq2'] = valueRow
    ## fisher test
    norCountsFile = os.path.join(outputDir, 'm6a_peak_diff.normalized.counts.txt')
    norCountsDict = defaultdict(dict)
    with open(norCountsFile, 'r') as f:
        sampleHeaderIndexDict = defaultdict(dict)
        hrow = f.readline().strip().split('\t')
        for i in range(len(hrow)):
            ## start from sample
            sampleHeaderIndexDict[i+1] = hrow[i].replace('X', '')
        for line in f:
            #start from peakid
            row = line.strip().split()
            peakid = row[0]
            if peakid not in peakDict:
                continue
            for i in range(1, len(row)):
                uniqBamIndex = sampleHeaderIndexDict[i]
                ipType = bamDict['indexToInfo'][uniqBamIndex]['ipType']
                treatType = bamDict['indexToInfo'][uniqBamIndex]['treatType']
                if treatType not in norCountsDict[peakid]:
                    norCountsDict[peakid][treatType] = defaultdict(list)
                norCountsDict[peakid][treatType][ipType].append(float(row[i]))
    #delete files from DESeq2Gene.R
    os.remove(os.path.join(outputDir, 'm6a_peak_diff.DESeq2.txt'))
    os.remove(os.path.join(outputDir, 'm6a_peak_diff.sig.DESeq2.txt'))
    os.remove(os.path.join(outputDir, 'm6a_peak_diff.normalized.counts.txt'))
    os.remove(os.path.join(outputDir, 'm6a_peak_diff.MA.pdf'))
    os.remove(os.path.join(outputDir, 'm6a_peak_diff.pvalue.histogram.pdf'))
    os.remove(os.path.join(outputDir, 'm6a_peak_diff.pvalueNorCounts.bar.pdf'))
    # run fisher exact test parallel
    countsPeakidList = sorted(norCountsDict.keys())
    fisherTestResList = []
    with Pool(processes=thread) as pool:
        for i in pool.imap_unordered(RunFisherExactTest, countsPeakidList):
            fisherTestResList.append(i)
    ## get fisher p-value and BH FDR
    fisherPeakidList = []
    fisherPvalList = []
    for result in fisherTestResList:
        peakid, pvalue = result
        fisherPeakidList.append(peakid)
        fisherPvalList.append(float(pvalue))
    padjList = rstats.p_adjust(FloatVector(fisherPvalList), method = 'BH')
    ## store pvalue and fdr
    fisherPvalDict = {}
    for i in range(len(fisherPeakidList)):
        peakid = fisherPeakidList[i]
        pvalue = fisherPvalList[i]
        padj = padjList[i]
        fisherPvalDict[peakid] = [pvalue, padj]


with open(args.output, 'w') as out:
    # delete fdr
    del nameRow[-3]
    if bool(args.bamcip):
        nameRow.extend(['mean.normalized.counts.DEseq2', 'diff.log2fc.DESeq2', 'diff.pvalue.fisher', 'diff.fdr.fisher'])
    if bool(args.expMtx):
        nameRow.extend(['gene.aveExp.control', 'gene.aveExp.treated'])
    if bool(args.degMtx):
        nameRow.extend(['gene.DE.log2fc', 'gene.DE.pval'])
    out.write('\t'.join(nameRow) + '\n')
    for peakid in peakidList:
        row = peakDict[peakid]['row']
        row[3] = row[3].split('|')[0]
        ## convert log.pval to pval
        row[-3] = 10 ** float(row[-3])
        ## delete fdr
        del row[-1]
        ## label the source of peak
        peakSource = peakDict[peakid]['source']
        row.append(peakSource)
        if bool(args.bamcip):
            if 'peakDiffDeseq2' not in peakDict[peakid]:
                continue
            if peakid not in fisherPvalDict:
                continue
            row += peakDict[peakid]['peakDiffDeseq2'] + fisherPvalDict[peakid]
        if args.expMtx:
            row += peakDict[peakid]['geneExp']
        if args.degMtx:
            row += peakDict[peakid]['geneDe']
        out.write('\t'.join(map(str, row)) + '\n')
