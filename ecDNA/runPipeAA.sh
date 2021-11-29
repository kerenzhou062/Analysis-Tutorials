#!/bin/bash
#SBATCH --job-name=ecdnaPrepareAA    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 1                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=60G                      # Amount of memory in GB
#SBATCH --time=300:10:00               # Time limit hrs:min:sec
#SBATCH --output=ecdnaPrepareAA.log               # Time limit hrs:min:sec

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.

function showHelp {
  echo -ne "usage: sbatch ecdnaPrepareAA.sh -n 1 -o <log> --mem <60G> "
  echo -ne "ecdnaPrepareAA.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -a | --anno: the annotations directory of ReadDepth (required by --call-readDepth) <str>
    -b | --bam: coordinate sorted bam (not with --read1 and --read2) <str>
    -c | --cnv: custom <str>
    -m | --memory: #GB of memory [60] <int>
    -o | --output: the output folder [./] <str>
    -p | --prefix: the output prefix [ecDNA] <str>
    -t | --thread: # of cpus <int>
    -r | --ref: reference name ['hg19','GRCh37','GRCh38'(default)] <str>
    -c | --chrom ) CHROM_FILE="$2"; shift 2 ;;
    --cngain: cutoff for CN gain considered by AA [5] <float>
    --cnsize: cutoff for CN interval size considered by AA [50000] <int>
    --downsample: cutoff for bam coverage downsampling during AA [10] <float>
    --read1: the input fastq-r1 file <str>
    --read2: the input fastq-r2 file <str>
    --call-readDepth: use cnvCallReadDepth.sh to call cnv instead of cnvkit <bool>
    --set-empty: delete whole folder before running <bool>
    --skip-cnv: skip cnv calling step <bool>
    --skip-mapping: skip reads mapping step <bool>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o ha:b:c:m:o:p:t:r:, --long anno:,bam:,chrom:,cnv:,memory:,output:,prefix:,thread:,ref:, \
  --long cngain:,cnsize:,conda:,downsample:,read1:,read2:, \
  --long help,call-readDepth,set-empty,skip-cnv,skip-mapping, \
  -- "$@"`


if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREAD=10
MEMORY=60
ANNO_DIR=
BAM_FILE=
PREFIX="ecDNA"
OUTPUT_DIR="./"
CNV_BED=""
REF_AA="GRCh38"
CHROM_FILE=
CN_GAIN=5
CN_SIZE=50000
DOWAN_SAMPLE=10
READ1=
READ2=
CALL_READ_DEPTH=false
SET_EMPTY=false
SKIP_CNV=false
SKIP_MAPPING=false

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -a | --anno ) ANNO_DIR="$2"; shift 2 ;;
    -b | --bam ) BAM_FILE="$2"; shift 2 ;;
    -c | --cnv ) CNV_BED="$2"; shift 2 ;;
    -m | --memory ) MEMORY="$2"; shift 2 ;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2 ;;
    -p | --prefix ) PREFIX="$2"; shift 2 ;;
    -r | --ref ) REF_AA="$2"; shift 2 ;;
    -t | --thread ) THREAD="$2"; shift 2 ;;
    --chrom ) CHROM_FILE="$2"; shift 2 ;;
    --cngain ) CN_GAIN="$2"; shift 2 ;;
    --cnsize ) CN_SIZE="$2"; shift 2 ;;
    --downsample ) DOWAN_SAMPLE="$2"; shift 2 ;;
    --read1 ) READ1="$2"; shift 2 ;;
    --read2 ) READ2="$2"; shift 2 ;;
    --call-readDepth ) CALL_READ_DEPTH=true; shift ;;
    --set-empty ) SET_EMPTY=true; shift ;;
    --skip-cnv ) SKIP_CNV=true; shift ;;
    --skip-mapping ) SKIP_MAPPING=true; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

## getopt end

echo "Running pipleline with following parameters:"
echo "ANNO_DIR=$ANNO_DIR"
echo "BAM_FILE=$BAM_FILE"
echo "THREAD=$THREAD"
echo "MEMORY=$MEMORY"
echo "PREFIX=$PREFIX"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo "CNV_BED=$CNV_BED"
echo "REF_AA=$REF_AA"
echo "CHROM_FILE=$CHROM_FILE"
echo "CN_GAIN=$CN_GAIN"
echo "CN_SIZE=$CN_SIZE"
echo "DOWAN_SAMPLE=$DOWAN_SAMPLE"
echo "READ1=$READ1"
echo "READ2=$READ2"
echo "CALL_READ_DEPTH=$CALL_READ_DEPTH"
echo "SET_EMPTY=$SET_EMPTY"
echo "SKIP_MAPPING=$SKIP_MAPPING"
echo "SKIP_CNV=$SKIP_CNV"
echo ""

## check required arguments
if [ -z $READ1 ] && [ ! -f $READ1 ] ; then
  if [ -z $BAM_FILE ]; then
    echo "--read1: Wrong! read1 NOT found!"
    exit 2
  fi
fi

if [ -z $READ2 ] && [ ! -f $READ2 ] ; then
  if [ -z $BAM_FILE ]; then
    echo "--read2: Wrong! read2 NOT found!"
    exit 2
  fi
fi

if [[ "$CALL_READ_DEPTH" == true ]]; then
  if [[ ! -f $CHROM_FILE ]]; then
    echo "--chrom: Wrong! chrom file NOT found!"
    exit 2
  fi
fi


# make directories
if [ "$SET_EMPTY" = true ]; then
  if [ "$SKIP_MAPPING" = false ]; then
    if [ ! -d $OUTPUT_DIR ]; then
      mkdir -p $OUTPUT_DIR
    else
      find $OUTPUT_DIR -type f | grep -Pv 'bam|.stderr' | xargs -I {} rm -rf {}
    fi
  fi
else
  if [ ! -d $OUTPUT_DIR ]; then
    mkdir -p $OUTPUT_DIR
  fi
fi

cd $OUTPUT_DIR

echo "Loading modules..."
module load AmpliconArchitect R/3.5.1

SKIP_MAP_ARG=""
if [ "$SKIP_MAPPING" = true ]; then
  SKIP_MAP_ARG="--skip_mapping"
fi

SORTED_BAM=false
if [ ! -z $READ1 ] && [ ! -z $READ2 ]; then
  SORTED_BAM="${OUTPUT_DIR}/${PREFIX}.cs.rmdup.bam"
  INPUT_ARG=" --fastqs $READ1 $READ2 --cngain $CN_GAIN --cnsize_min $CN_SIZE "
else
  if [ ! -z $BAM_FILE ] && [ -f $BAM_FILE ]; then
    echo "-b|--bam: Wrong! bam file NOT found!"
    exit 2
  else
    INPUT_ARG=" --sorted_bam $BAM_FILE --cngain $CN_GAIN --cnsize_min $CN_SIZE "
    SORTED_BAM=$BAM_FILE
  fi
fi

if test -f "$SORTED_BAM"; then
  SKIP_MAPPING=true
  INPUT_ARG=" --sorted_bam $SORTED_BAM --cngain $CN_GAIN --cnsize_min $CN_SIZE "
  echo "$SORTED_BAM file found! Automatically skip genome mapping!"
fi

CNV_ARG=""

if [ "$CNV_BED" != "" ] && [ -f "$CNV_BED" ]; then
  CNV_ARG=" --cnv_bed ${CNV_BED} "
fi

if [ "${CALL_READ_DEPTH}" = true ]; then
  if [ ! -z $ANNO_DIR ] && [ -d $ANNO_DIR ]; then
    if [ "$SKIP_MAPPING" = false ]; then
      echo "Mapping reads to genome by bwa:${REF_AA}"
      runBwa.py -s $PREFIX --thread $THREAD --memory ${MEMORY} --ref $REF_AA \
        --fastqs $READ1 $READ2 -o $OUTPUT_DIR > runBwa.log 2>&1
    fi

    if [ "$SORTED_BAM" != false ]; then
      READ_DEPTH_RES="${OUTPUT_DIR}/readDepth_output"
      CNV_BED="${READ_DEPTH_RES}/${PREFIX}.CNVs_GAIN.bed"
      if [ "$SKIP_CNV" = false ]; then
        echo "Running cnvCallReadDepth.sh to generate CNVs BED"
        ## run cnvCallReadDepth.sh to generate CNVs BED
        cnvCallReadDepth.sh --thread $THREAD --anno $ANNO_DIR --bam $SORTED_BAM -o $READ_DEPTH_RES -p $PREFIX \
          --cngain $CN_GAIN --cnsize $CN_SIZE --chrom ${CHROM_FILE} --paired > ${OUTPUT_DIR}/${PREFIX}.readDepth.log 2>&1
        CNV_ARG=" --cnv_bed ${CNV_BED} "
      else
        echo "Skip cnv calling with cnvCallReadDepth.sh"
        CNV_ARG=" --skip_cnv "
      fi

      ## reset INPUT_ARG
      INPUT_ARG=" --sorted_bam $SORTED_BAM "
    fi
  else
    echo "-a|--anno: Wrong! Not a folder!"
    exit 2
  fi
else
  if [ "$SKIP_MAPPING" = true ]; then
    SKIP_MAP_ARG="--skip_mapping"
  fi
  if [ "$SKIP_CNV" = true ]; then
    CNV_ARG=" --skip_cnv "
  fi
fi

echo "Running pipeAA.py..."
echo "INPUT_ARG=${INPUT_ARG}"
echo "CNV_ARG=${CNV_ARG}"
echo "SKIP_MAP_ARG=${SKIP_MAP_ARG}"

pipeAA.py -s $PREFIX  --thread $THREAD --memory ${MEMORY} ${INPUT_ARG} ${CNV_ARG} --run_AA \
  --downsample $DOWAN_SAMPLE --ref $REF_AA ${SKIP_MAP_ARG} -o $OUTPUT_DIR 

echo "All jobs done!"
