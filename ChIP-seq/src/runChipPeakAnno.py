#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict
import shutil
from multiprocessing import Pool
import subprocess
from PubAlbum import Anno

# usage: runChipPeakAnno.py or runChipPeakAnno.py <peak dir> -genome hg38

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='threads used for homer')
parser.add_argument('-fdr', action='store', type=str,
                    default="0.05,0.05", help='q-value cutoff(narrow,broad)')
parser.add_argument('-fold', action='store', type=str,
                    default="2,2", help='fold_enrichment cutoff(narrow,broad)')
parser.add_argument('-genome', action='store', type=str,
                    default='hg38', help='genome version (homer)')
parser.add_argument('-gtf', action='store', type=str,
                    default='hg38v32',
                    help='gtf annotation build or file("None" for homer default)')
parser.add_argument('-input', action='store', type=str,
                    default='./peak', help='The main peak folder')
parser.add_argument('-matrix', action='store', type=str,
                    help='The input matrix (when sbatch activate)')
parser.add_argument('-memory', action='store', type=str,
                    default='100G', help='memory used for sbatch')
parser.add_argument('-output', action='store', type=str,
                    default='./peakAnnoRes', help='annoation output directory')
parser.add_argument('-part', action='store', type=str,
                    default='all', help='partition of slurm server')
parser.add_argument('-pval', action='store', type=str,
                    default="1e-4,1e-2", help='p-value cutoff(narrow,broad)')
parser.add_argument('-sbatch', action='store', 
                    default=False, help='activate sbatch mode')
parser.add_argument('-size', action='store_true', type=str,
                    default='8,10,12', help='-len parameter (homer)')
parser.add_argument('-top', action='store', type=int,
                    default=3,
                    help='Top motif for annotatePeaks.pl')
parser.add_argument('-txdb', action='store', type=str,
                    default='hg38v32',
                    help='TxDb object build or file (.txdb)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    print("Running runChipPeakAnno.py with defaultdict parameters...")


def runHomer(matrixDict, peakFile, genome, gtf, top, size, cpu):
    # annotate genes with annotatePeaks.pl
    def anno(peak, output):
        # find motifs with findMotifsGenome.pl
        os.chdir(output)
        tmpDir = os.path.join(output, 'tmp')
        motifDir = ''
        motfiCommand = 'findMotifsGenome.pl {peak} {genome} {} -p {cpu} '
        # annotate peaks with annotatePeaks.pl
        gtfPram = '-gtf {gtf}'.format(**vars())
        if gtf == 'None':
            gtfPram = ''
        annoPeakResFile = 'homer_' + prefix + '.anno'
        annoCommand = 'annotatePeaks.pl {peakFile} {genome} \
            {gtfPram} > {annoPeakResFile}'.format(**vars())

    basename = os.path.basename(peakFile)
    prefix = os.path.splitext(basename)[0]
    outputDir = os.path.splitext(peakFile)[0]



# public arguments
threadNum = args.cpu
genome = args.genome
gtf = args.gtf
memory = args.memory
partition = args.part
sbatchFlag = args.sbatch
top = args.top
folds = list(map(float, args.fold.split(',')))
pvals = list(map(float, args.pval.split(',')))
fdrs = list(map(float, args.fdr.split(',')))
narrowPeak = 'narrowPeak'
broadPeak = 'broadPeak'
anno = Anno()
gtf = anno.gtf(args.gtf)
basepath = os.path.realpath(args.input)
output = os.path.realpath(args.output)

# based arguments
bashDir = os.path.join(basepath, 'runBash')
logDir = os.path.join(bashDir, 'log')
matrixFile = os.path.join(output, 'peakFile.matrix')

# ================================
# sbatch
# ================================

sbatchAnnoTemplate = '''#!/bin/bash
#SBATCH --job-name=ChipPeakAnno_{baseExpName}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kzhou@coh.org
#SBATCH -n {threadNum}
#SBATCH -N 1-1
#SBATCH -p all
#SBATCH --mem={memory}
#SBATCH --time=72:00:00
#SBATCH --output={logPath}

# ================================
# load required module
# ================================
echo -e "load modules..."
module load bedtools/2.29.0
module load Python/3.6.5
module load Anaconda2/4.2.0

# ===========step 0a==============
# public bariables
# ================================
echo -e "\\nPreparing basic variables..."

BASE="{basepath}"
THREADS={threadNum}
GENOME="{genome}"
MATRIX="{matrix}"
GTF="{gtf}"

# =========peak annotation with homer ===========
# Parallel running
# =======================
echo "Start peak annotation with homer..."

runChipPeakAnno.py -cpu ${{THREADS}} \\
  -matrix ${{MATRIX}} \\
  -genome ${{GENOME}} \\
  -gtf ${{GTF}} \\
  -pval {pval} \\
  -output ${{OUTPUT}} \\
  -sbatch ${{RANK}}

echo "IDR done..."
'''

if sbatchFlag is False:
    # try to delete dirs
    try:
        shutil.rmtree(output)
    except FileNotFoundError:
        pass
    os.makedirs(output, exist_ok=True)
    # file pattern: Kas1_control_NA_rep1_final_peaks.xls,
    # Kas1_control_NA_rep1_final.narrowPeak,
    # Kas1_control_NA_rep1_final.broadPeak
    extsList = [narrowPeak, broadPeak]
    peakFileList = sorted([f for ext in extsList for f in glob(
        os.path.join(basepath, '**', '*.' + ext), recursive=True)])
    cleanCompile = re.compile(r'(?:_peaks|_with\w+)\..+$')
    repBoolDict = defaultdict(dict)
    for peakFile in peakFileList:
        # judgement of replicates
        basename = os.path.basename(peakFile)
        # expName: Kas1_control_NA_pooled
        expName = cleanCompile.sub('', basename)
        if re.search(r'/idr/', basename, re.IGNORECASE):
            continue
        if re.search(r'pr', basename, re.IGNORECASE):
            continue
        if expName in repBoolDict:
            continue
        if bool(re.search(r'\/narrow\/', peakFile)):
            # prefixDir: path/narrow/
            prefixDir = os.path.split(os.path.split(peakFile)[0])[0]
            dirList = [name for name in os.listdir(
                prefixDir) if os.path.isdir(os.path.join(prefixDir, name))]
            rep2List = list(filter(lambda x: re.search(r'rep2', x), dirList))
            if bool(rep2List):
                repBoolDict[expName] = True
            else:
                repBoolDict[expName] = False
    
    matrixDict = defaultdict(dict)
    for peakFile in peakFileList:
        # get basename and information
        basename = os.path.basename(peakFile)
        expName = cleanCompile.sub('', basename)
        # remove non-necessary peak
        if re.search(r'/idr/', basename, re.IGNORECASE):
            continue
        if bool(re.search(r'pr', basename, re.IGNORECASE)):
            continue
        if expName not in repBoolDict:
            continue
        else:
            if repBoolDict[expName] is False:
                if bool(re.search(r'pool', basename, re.IGNORECASE)):
                    continue
        ext = os.path.splitext(basename)[-1]
        # label peak file
        if bool(re.search(r'IDR', basename, re.IGNORECASE)):
            idrLabel = True
        else:
            idrLabel = False
        if bool(re.search(rf'{narrowPeak}', ext)):
            peakLabel = narrowPeak
        elif bool(re.search(rf'{broadPeak}', ext)):
            peakLabel = broadPeak
        else:
            continue
        matrixDict[peakFile]['peakLabel'] = peakLabel
        matrixDict[peakFile]['idrLabel'] = idrLabel
        matrixDict[peakFile]['prefix'] = os.path.splitext(basename)[0]
        outputDir = os.path.join(output, peakLabel)
        if idrLabel:
            outputDir = os.path.join(output, 'idr_' + peakLabel)
        matrixDict[peakFile]['output'] = outputDir
        matrixDict[peakFile]['dest'] = os.path.join(outputDir, basename)
        
    with open(matrixFile, 'w') as matrixO:
        # source peakFile, destination, narrow, True, 
        # HepG2_CTCF_NA_pooled_withRep.IDR.filt, output dir
        row = ['source', 'dest', 'peakLabel', 'idrLabel', 'prefix', 'output']
        matrixO.write('\t'.join(row) + '\n')
        for peakFile in sorted(matrixDict.keys()):
            dest = matrixDict[peakFile]['dest']
            peakLabel = matrixDict[peakFile]['peakLabel']
            idrLabel = matrixDict[peakFile]['idrLabel']
            prefix = matrixDict[peakFile]['prefix']
            output = matrixDict[peakFile]['output']
            idrLabel = 1 if idrLabel else 0
            row = [peakFile, dest, peakLabel, idrLabel, prefix, output]
            matrixO.write('\t'.join(row) + '\n')
    # generate sbatch scripts and runBash script
    runSbatchScript = os.path.join(bashDir, "runChipPeakAnno.sbatch")
    sbatchScript = os.path.join(bashDir, 'runPeaksAnnotation.sh')
    logPath = os.path.join(logDir, 'runPeaksAnnotation.log')
    with open(runSbatchScript, 'w') as sbatchO, open(sbatchScript, 'w') as sbatchScriptO:
        sbatchO.write('#!/bin/sh\n')
        sbatchO.write('BASE={0}\n\n'.format(bashDir))
        sbatchO.write('#sbatch for annotation\n')
        sbatchO.write('sbatch $BASE/{sbatchScript}\n'.format(**vars()))
        sbatchCont = sbatchAnnoTemplate.format(**vars())
        sbatchScriptO.write(sbatchCont)
    os.chmod(sbatchScript, 0o755)
    os.chmod(runSbatchScript, 0o755)
else:
    matrixDict = defaultdict(dict)
    # run in sbatch mode
    with open(matrix, 'r') as f:
        __ = f.readline()
        for line in f:
            row = line.strip().split('\t')
            source = row[0]
            dest = row[1]
            peakLabel = row[2]
            idrLabel = int(row[3])
            prefix = row[4]
            output = row[5]
            matrixDict[dest]['peakLabel'] = peakLabel
            matrixDict[dest]['idrLabel'] = True if idrLabel else False
            matrixDict[dest]['prefix'] = prefix
            matrixDict[dest]['output'] = output
            # copy peak files
            os.makedirs(output, exist_ok=True)
            cpfile = shutil.copyfile(source, dest)
    # annotate peak in parallel
    sys.stderr.write("Running Homer...")
    pool = Pool(processes=args.cpu)
    for peakFile in sorted(matrixDict.keys()):
        peakLabel = matrixDict[peakFile]['peakLabel']
        idrLabel = matrixDict[peakFile]['idrLabel']
        prefix = matrixDict[peakFile]['prefix']
        output = matrixDict[peakFile]['output']
        # try to make dirs
        try:
            shutil.rmtree(output)
        except FileNotFoundError:
            pass
        os.makedirs(output, exist_ok=True)
        # runMacs2(ip, control, name, other,
        #    output, extsize, gsize, shift, pval, dFormat, dtype, rank, chrsize)
        result = pool.apply_async(runHomer, args=(matrixDict, matrixFile))
    pool.close()
    pool.join()

    endtime = datetime.datetime.now()
    collapsed = (endtime - starttime).seconds
    sys.stderr.write("All jobs done (Homer)!")
    sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
