#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict
from itertools import combinations
from PubAlbum import Anno

#usage: runExomePeakBash.py or runExomePeakBash.py <bam dir>

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-control', action='store', type=str,
                    required=True,
                    default='shCont', help='keyword for control samples')
parser.add_argument('-gtf', action='store', type=str,
                    required=True,
                    help='gtf annotation file')
parser.add_argument('-grepKept', action='store', type=str,
                    help='regex for keeping files')
parser.add_argument('-grepExpel', action='store', type=str,
                    help='regex for filtering files')
parser.add_argument('-input', action='store', type=str,
                    default='./',
                    help='input bam directory')
parser.add_argument('-log', action='store', type=str,
                    default='./',
                    help='log directory')
parser.add_argument('-memory', action='store', type=str,
                    default='40G',
                    help='memory used for sbatch')
parser.add_argument('-output', action='store', type=str,
                    default='./',
                    help='output folder')
parser.add_argument('-package', action='store', type=str,
                    choices=['exomePeak', 'QNB'],
                    default='QNB',
                    help='package for detection of differential methylated peaks')
parser.add_argument('-prefix', action='store', type=str,
                    default='result', help='prefix for output')
parser.add_argument('-script', action='store', type=str,
                    default='./',
                    help='directory for storing scripts')
parser.add_argument('--rdata', action='store_true',
                    default=False,
                    help='load exomePeakRes.RData from output directory instead')
parser.add_argument('-treat', action='store', type=str,
                    required=True,
                    help='keyword for treatment samples')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

##public arguments
memory = args.memory
gtf = args.gtf
log = args.log
package = args.package
prefix = args.prefix
output = args.output
script = args.script
basepath = os.path.realpath(args.input)

if bool(args.grepKept):
    kept = True
    regex = re.compile(r'{0}'.format(args.grepKept))
elif bool(args.grepExpel):
    kept = False
    regex = re.compile(r'{0}'.format(args.grepExpel))
else:
    regex = False

bamFiles = sorted(glob(os.path.join(basepath, '**', '*.bam'), recursive=True))

controlBamDict = defaultdict(list)
treatBamDict = defaultdict(list)
cRegex = re.compile(r'{0}'.format(args.control))
tRegex = re.compile(r'{0}'.format(args.treat))

for bamFile in bamFiles:
    fileName = os.path.split(bamFile)[-1]
    if bool(regex):
        if kept:
            if bool(regex.search(fileName)) is False:
                continue
        else:
            if bool(regex.search(fileName)) is True:
                continue
    if bool(cRegex.search(fileName)) is True:
        if re.search(r'IP', fileName, flags=re.IGNORECASE):
            controlBamDict['IP'].append(bamFile)
        else:
            controlBamDict['input'].append(bamFile)
    if bool(tRegex.search(fileName)) is True:
        if re.search(r'IP', fileName, flags=re.IGNORECASE):
            treatBamDict['IP'].append(bamFile)
        else:
            treatBamDict['input'].append(bamFile)

if args.rdata:
    rdata = 'TRUE'
else:
    rdata = 'FALSE'

qnbRTemplate = '''#!/usr/bin/env Rscript
#QNB Script
#R script
#Define parameters and load library
library("exomePeak")
library("QNB")
library("GenomicFeatures")
library("Rsamtools")

getCounts <- function(bamList, peak, exonRangesObj, countsObj) {{
  require("exomePeak")
  require("QNB")
  require("GenomicFeatures")
  require("GenomicRanges")
  require("Rsamtools")
  for (i in 1:length(bamList)) {{
    bam <- bamList[i]
    aligns <- readGAlignments(bam)
    id <- countOverlaps(aligns,exonRangesObj)
    txFilteredAligns <- aligns[id>0]
    countsObj[,i] <- countOverlaps(peak, txFilteredAligns)
  }}
  return(countsObj)
}}


GENE_ANNO_GTF = "{gtf}"

resDir = file.path('{output}', '{prefix}')
dir.create(resDir, showWarnings = FALSE)
setwd(resDir)

UNTREATED_IP_BAM <- c({controlIpBams})
UNTREATED_INPUT_BAM <- c({controlInputBams})
TREATED_IP_BAM <- c({treadIpBams})
TREATED_INPUT_BAM <- c({treadInputBams})

IP_BAM <- c({controlIpBams}, {treadIpBams})
INPUT_BAM <- c({controlInputBams}, {treadInputBams})

# running exomePeak
rdataFlag <- {rdata}
if (isTRUE(rdataFlag)) {{
  load("exomePeakRes.RData")
}}else{{
  res <- exomepeak(GENE_ANNO_GTF=GENE_ANNO_GTF, IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM)
  save(res, file='exomePeakRes.RData')
}}

#get reads count
peak <- res$all_peaks
untreated_ip <- matrix(0,nrow=length(peak),ncol=length(UNTREATED_IP_BAM))
untreated_input <- matrix(0,nrow=length(peak),ncol=length(UNTREATED_INPUT_BAM))
treated_ip <- matrix(0,nrow=length(peak),ncol=length(TREATED_IP_BAM))
treated_input <- matrix(0,nrow=length(peak),ncol=length(TREATED_INPUT_BAM))

txdb <- makeTxDbFromGFF(file=GENE_ANNO_GTF, format="gtf")
exonRanges <- exonsBy(txdb, "tx")

#get ip reads count
#control_IP
untreated_ip <- getCounts(UNTREATED_IP_BAM, peak, exonRanges, untreated_ip)

#treadIP
treated_ip <- getCounts(TREATED_IP_BAM, peak, exonRanges, treated_ip)

#get input reads count
#control_Input
untreated_input <- getCounts(UNTREATED_INPUT_BAM, peak, exonRanges, untreated_input)

#treadInput
treated_input <- getCounts(TREATED_INPUT_BAM, peak, exonRanges, treated_input)

#differential RNA methylation analysis
result = qnbtest(untreated_ip, treated_ip, untreated_input, treated_input)

##paste files
setwd(resDir)

allPeakFile <- "exomePeak_output/peak.xls"
qnbResult <- "dif_meth.xls"
awkMain <- '{{if(FNR==1){{$0=gensub(/"/,"","g",$0)}}print}}'
system(paste("paste", "-d", "$'\\t'", allPeakFile, qnbResult, "|", "awk", "'", awkMain, "'", ">", "{prefix}.diff.xls", sep=" "))

sessionInfo()
'''


exompeakRTemplate = '''#!/usr/bin/env Rscript
#exomepeak Script
#R script
#Define parameters and load library
library("exomePeak")

base = "{output}"
setwd(base)

GENE_ANNO_GTF = "{gtf}"
# differential peak calling

IP_BAM=c({controlIpBams})
INPUT_BAM=c({controlInputBams})
TREATED_IP_BAM=c({treadIpBams})
TREATED_INPUT_BAM=c({treadInputBams})

resDir = file.path('{output}', '{prefix}')
dir.create(resDir, showWarnings = FALSE)
setwd(resDir)

result = exomepeak(GENE_ANNO_GTF=GENE_ANNO_GTF, 
    IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM, 
    TREATED_IP_BAM=TREATED_IP_BAM, TREATED_INPUT_BAM=TREATED_INPUT_BAM)

##rename and remove files
file.rename("./exomePeak_output/con_sig_peak.bed", "./{prefix}.con_sig_peak.bed")
file.rename("./exomePeak_output/con_sig_peak.xls", "./{prefix}.con_sig_peak.xls")
file.rename("./exomePeak_output/peak.bed", "./{prefix}.peak.bed")
file.rename("./exomePeak_output/peak.xls", "./{prefix}.peak.xls")
file.rename("./exomePeak_output/sig_peak.bed", "./{prefix}.sig_peak.bed")
file.rename("./exomePeak_output/sig_peak.xls", "./{prefix}.sig_peak.xls")
file.rename("./exomePeak_output/exomePeak.Rdata", "./{prefix}.exomePeak.Rdata")
file.remove("./exomePeak_output")

sessionInfo()
'''

sbatchTemplate = '''#!/bin/sh
#SBATCH --job-name={rScriptFile}    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail  
#SBATCH -n 1                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem={memory}                      # Amount of memory in GB
#SBATCH --time=100:00:00               # Time limit hrs:min:sec
#SBATCH --output={log}/{rScriptFile}.log   # Standard output and error log

BASE={script}

Rscript $BASE/{rScriptFile}
'''

controlIpBams = ', '.join(map(lambda x:'"{}"'.format(x), sorted(controlBamDict['IP'])))
controlInputBams = ', '.join(map(lambda x:'"{}"'.format(x), sorted(controlBamDict['input'])))
treadIpBams = ', '.join(map(lambda x:'"{}"'.format(x), sorted(treatBamDict['IP'])))
treadInputBams = ', '.join(map(lambda x:'"{}"'.format(x), sorted(treatBamDict['input'])))

rScriptFile = "{prefix}.{package}.DifferPeak.r".format(**vars())
runRScript = os.path.join(script, rScriptFile)

with open(runRScript, 'w') as runRScriptO:
    if package == 'QNB':
        runRScriptO.write(qnbRTemplate.format(**vars()))
    else:
        runRScriptO.write(exompeakRTemplate.format(**vars()))
os.chmod(runRScript, 0o755)

runBashScript = os.path.join(script, prefix + ".DifferPeak.sh")
with open(runBashScript, 'w') as runBashScriptO:
    runBashScriptO.write(sbatchTemplate.format(**vars()))

os.chmod(runBashScript, 0o755)
