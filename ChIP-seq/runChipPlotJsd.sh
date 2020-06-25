#!/bin/bash

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
# ref:https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c

function showHelp {
  echo -ne "usage: runChipPlotJsd.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -i | --input: input directory (eg. alignment/HepG2_CTCF_WT, contains IP&input folder) <str>
    -p | --prefix: the prefix of output JSD.PNG file (eg. HepG2_CTCF_WT) <str>
    -q | --quality: the mapping quality cutoff [30] <str>
    -t | --thread: # of cpus [10] <int>
    --load: load necessary modules on slurm server <int>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o h:i:p:q:t:, --long help,thread:,input:,prefix:,quality:, \
  --long load, \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREADS=10
MAPQ_THRESH=30
PREFIX=
INPUT_DIR=
LOAD_FLAG=false

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -i | --input ) INPUT_DIR="$2"; shift 2 ;;
    -p | --prefix ) PREFIX="$2"; shift 2 ;;
    -q | --quality ) MAPQ_THRESH="$2"; shift 2 ;;
    -t | --thread ) THREADS="$2"; shift 2 ;;
    --load ) LOAD_FLAG=true; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done
## getopt end

## check required arguments
if [ -z "$INPUT_DIR" ] || [ ! -d "$INPUT_DIR" ]; then
  echo "--input: Wrong! input directory NOT found!"
  exit 2
fi

if [ -z "$PREFIX" ]; then
  echo "--prefix: Wrong! prefix for output not SET!"
  exit 2
fi

echo "Running pipleline with following parameters:"
echo "THREADS=$THREADS"
echo "INPUT_DIR=$INPUT_DIR"
echo "MAPQ_THRESH=$MAPQ_THRESH"
echo "PREFIX=$PREFIX"
echo ""

if [ "$LOAD_FLAG" = true ]; then
  echo -e "load modules..."
  module load picard/2.21.1
  module load bedtools/2.29.0
  module load Python/3.6.5
  module load Anaconda2/4.2.0
fi

echo "Change directory to ${INPUT_DIR}"
cd "${INPUT_DIR}"

FILE_PREFIX="${INPUT_DIR}/${PREFIX}"

# ========step 2f============
# Calculate Jensen-Shannon distance (JSD)
# =======================
# output:
# ● JSD Plot ${JSD_PLOT}
# ● JSD log ${JSD_LOG}
echo "Step 2f starting..."

JSD_LOG="${FILE_PREFIX}.JSD.log"
JSD_PLOT="${FILE_PREFIX}.JSD.png"

# BAMs are blacklist-filtered first for each replicate and control
# FINAL_BAM_PREFIX="${FILE_PREFIX}.filt.srt.nodup"
# FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored

BAMS=""
LABELS=""
for bam in `find ./ -maxdepth 2 -type f -name "*filt.srt.nodup.bfilt.bam"| grep -v 'pooled'| sort`;
do
  BAMS=`echo "${BAMS}"" ""${bam}"`
  BASE_NAME=$(basename ${bam%%.filt*.bam})
  LABELS=`echo "${LABELS}"" ""${BASE_NAME}"`
done

echo "Plotting Jensen-Shannon distance"

C_ALL=en_US.UTF-8 LANG=en_US.UTF-8 plotFingerprint -b ${BAMS} \
  --labels ${LABELS} --outQualityMetrics ${JSD_LOG} \
  --minMappingQuality ${MAPQ_THRESH} -T "Fingerprints of different samples" \
  --numberOfProcessors ${THREADS} --plotFile ${JSD_PLOT}

echo "Delete *filt.srt.nodup.bfilt.bam!"
find "$INPUT_DIR" -maxdepth 2 -type f -name "*filt.srt.nodup.bfilt.bam" | xargs -I {} rm -rf {}
find "$INPUT_DIR" -maxdepth 2 -type f -name "*filt.srt.nodup.bfilt.bam.bai" | xargs -I {} rm -rf {}
#
echo "Step 2f done!"
