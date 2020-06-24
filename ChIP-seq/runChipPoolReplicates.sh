#!/bin/bash

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
# ref:https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c

function showHelp {
  echo -ne "usage: runChipPoolReplicates.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -b | --blacklist: input genome blacklist bed, which consists of centromere, repeat regions, etc. <str>
    -i | --input: input directory (eg. alignment/HepG2_CTCF_WT/IP) <str>
    -p | --prefix: the prefix of dataset (eg. HepG2_CTCF_WT_IP) <str>
    --load: load necessary modules on slurm server <int>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o hb:i:p:, --long help,input:,blacklist:,prefix:, \
  --long load, \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
BLACKLIST=
INPUT_DIR=
PREFIX=
LOAD_FLAG=false

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -b | --blacklist ) BLACKLIST="$2"; shift 2 ;;
    -i | --input ) INPUT_DIR="$2"; shift 2 ;;
    -p | --prefix ) PREFIX="$2"; shift 2 ;;
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
echo "BLACKLIST=$BLACKLIST"
echo "INPUT_DIR=$INPUT_DIR"
echo "PREFIX=$PREFIX"
echo ""

if [ "$LOAD_FLAG" = true ]; then
  echo -e "load modules..."
  module load picard/2.21.1
  module load bedtools/2.29.0
  module load Python/3.6.5
  module load Anaconda2/4.2.0
fi

# ===========step 0==============
# public bariables
# ================================
echo -e "\nRunning runChipPoolReplicates."
echo -e "\nPreparing basic variables..."

DATASET_PREFIX="${INPUT_DIR}/${PREFIX}"

echo "Change wd directory to ${INPUT_DIR}"
cd ${INPUT_DIR}

# =========step 2d===========
# Create pooled datasets
# =======================
# output:
# ● Pooled tagAlign file ${POOLED_TA_FILE}
# ● 2 pooled-pseudoreplicate tagAlign files
echo "Step 2d starting..."
echo "Pooling TA files of real replicates..."

POOLED_TA_FILE="${DATASET_PREFIX}.pooled.tagAlign.gz"
FINAL_TA_FILES=""
for i in `find ./ -maxdepth 1 -type f -links +1 -name "*.final.tagAlign.gz"|grep -v 'pooled'|sort`;
do
  FINAL_TA_FILES=`echo "${FINAL_TA_FILES}"" ""${i}"`
done

zcat -f ${FINAL_TA_FILES} | gzip -nc > ${POOLED_TA_FILE}

# ========================
# Create pooled pseudoreplicates
# =======================
# PR1_TA_FILE="${PR_PREFIX}.pr1.tagAlign.gz"
# PR2_TA_FILE="${PR_PREFIX}.pr2.tagAlign.gz"

echo "Pooling TA files of pseudo-replicates..."

PPR1_TA_FILE="${DATASET_PREFIX}.pooled.pr1.tagAlign.gz"
PPR2_TA_FILE="${DATASET_PREFIX}.pooled.pr2.tagAlign.gz"

PPR1_TA_FILES=""
for i in `find ./ -maxdepth 1 -type f -links +1 -name "*.pr1.tagAlign.gz"|grep -v 'pooled'|sort`;
do
  PPR1_TA_FILES=`echo "${PPR1_TA_FILES}"" ""${i}"`
done

PPR2_TA_FILES=""
for i in `find ./ -maxdepth 1 -type f -links +1 -name "*.pr2.tagAlign.gz"|grep -v 'pooled'|sort`;
do
  PPR2_TA_FILES=`echo "${PPR2_TA_FILES}"" ""${i}"`
done

zcat -f ${PPR1_TA_FILES} | gzip -nc > ${PPR1_TA_FILE}
zcat -f ${PPR2_TA_FILES} | gzip -nc > ${PPR2_TA_FILE}


echo "Step 2d done!"

# ========step 2f============
# Calculate Jensen-Shannon distance (JSD)
# =======================
# output:
# ● JSD Plot ${JSD_PLOT}
# ● JSD log ${JSD_LOG}
echo "Step 2f starting..."
# BAMs are blacklist-filtered first for each replicate and control
# FINAL_BAM_PREFIX="${DATASET_PREFIX}.filt.srt.nodup"
# FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored

echo "Filtering bam with blacklist bed file..."
NODUP_BFILT_BAMS=""
LABELS=""
for i in `find ./ -maxdepth 1 -type f -links +1 -name "*.filt.srt.nodup.bam"| grep -v 'pooled'| sort`;
do
  BTILT_PREFIX=${i%%.bam}
  NODUP_BFILT_BAM=${BTILT_PREFIX}.bfilt.bam
  FINAL_BAM_INDEX_FILE=${BTILT_PREFIX}.bfilt.bai
  bedtools intersect -nonamecheck -v -abam ${i} -b ${BLACKLIST} > ${NODUP_BFILT_BAM}
  samtools index ${NODUP_BFILT_BAM} ${FINAL_BAM_INDEX_FILE}
done

echo "Step 2f done!"

find "$INPUT_DIR" -maxdepth 1 -type f -links +1 -name "*.filt.*.bam" | grep -v '.filt.srt.nodup.bfilt' | xargs -I {} rm -rf {}
find "$INPUT_DIR" -maxdepth 2 -type f -name "*.filt.*.bam" | grep -v '.filt.srt.nodup.bfilt' | xargs -I {} rm -rf {}
find "$INPUT_DIR" -maxdepth 2 -type f -name "*.trim.*.bam" | grep -v '.filt.srt.nodup.bfilt' | xargs -I {} rm -rf {}
