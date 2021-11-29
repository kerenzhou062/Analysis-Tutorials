#!/bin/bash

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
# ref:https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c

function showHelp {
  echo -ne "usage: runChipPoolReplicates.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -b | --blacklist: input genome blacklist bed, which consists of centromere, repeat regions, etc. <str>
    -g | --gsize: the genome size <str>
    -o | --output: the output directory <str>
    -p | --prefix: the prefix of dataset (eg. HepG2_CTCF_WT) <str>
    -t | --type: peak type (narrowPeak|broadPeak) <str>
    --peak1: input peak1 <str>
    --peak2: input peak2 <str>
    --load: load necessary modules on slurm server <int>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o hb:g:p:o:t:, --long help,blacklist:,gsize:,output:,prefix:,\
  --long type:,peak1:,peak2:, \
  --long load, \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
PEAK_TYPE="narrowPeak"
LOAD_FLAG=false

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -b | --blacklist ) BLACKLIST="$2"; shift 2 ;;
    -g | --gsize ) GENOME_SIZE="$2"; shift 2 ;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2 ;;
    -p | --prefix ) PREFIX="$2"; shift 2 ;;
    -t | --type ) PEAK_TYPE="$2"; shift 2 ;;
    --peak1 ) PEAK1="$2"; shift 2 ;;
    --peak2 ) PEAK2="$2"; shift 2 ;;
    --load ) LOAD_FLAG=true; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done
## getopt end

## check required arguments
if [ -z "$OUTPUT_DIR" ]; then
  echo "--output: Wrong! output directory NOT found!"
  exit 2
fi

if [ -z "$PEAK_TYPE" ]; then
  echo "--type: Wrong! peak type not SET!"
  exit 2
fi

if [ "$PEAK_TYPE" != "narrowPeak" ] && [ "$PEAK_TYPE" != "broadPeak" ]; then
  echo "--type: Wrong! peak type not SET!"
  exit 2
fi

if [ -z "$PREFIX" ]; then
  echo "--prefix: Wrong! prefix for output not SET!"
  exit 2
fi

echo "Running pipleline with following parameters:"
echo "BLACKLIST=$BLACKLIST"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo "PREFIX=$PREFIX"
echo "PEAK1=$PEAK1"
echo "PEAK2=$PEAK2"
echo ""

if [ "$LOAD_FLAG" = true ]; then
  echo -e "load modules..."
  module load bedtools/2.29.0
fi

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

FILE_PREFIX="${OUTPUT_DIR}/${PREFIX}"

cd ${OUTPUT_DIR}

##overlap
echo "Start to overlap peaks..."

bedtools intersect -a $PEAK1 \
  -b $PEAK2 \
  -nonamecheck \
  > "${FILE_PREFIX}.overlap.${PEAK_TYPE}"

bedtools intersect -a $PEAK1 \
  -b $PEAK2 \
  -nonamecheck | \
  bedtools intersect -a stdin -b ${BLACKLIST} \
  -v -nonamecheck \
  > "${FILE_PREFIX}.overlap.filt.${PEAK_TYPE}"

if [ ! -z "$GENOME_SIZE" ] && [ -f "$GENOME_SIZE" ]; then
  ##fisher test
  echo "Start fisher test..."
  peak1Name=$(basename $PEAK1)
  peak2Name=$(basename $PEAK2)
  bedtools sort -i $PEAK1 -g $GENOME_SIZE > $peak1Name
  bedtools sort -i $PEAK2 -g $GENOME_SIZE > $peak2Name
  bedtools fisher -a $peak1Name -b $peak2Name -g $GENOME_SIZE \
    > "${FILE_PREFIX}.${PEAK_TYPE}.fisher.txt"
  rm -f $peak1Name $peak2Name
fi
