#!/bin/bash
#SBATCH --job-name=csv_call_by_readDepth    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 1                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=60G                      # Amount of memory in GB
#SBATCH --time=120:10:00               # Time limit hrs:min:sec
#SBATCH --output=csv_call_by_readDepth.log               # Time limit hrs:min:sec

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.

function showHelp {
  echo -ne "usage: sbatch cnvCallReadDepth.sh -n 1 -o <log> --mem <60G> "
  echo -ne "cnvCallReadDepth.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -a | --anno: the annotations directory of ReadDepth <str>
    -b | --bam: the input bam file <str>
    -c | --chrom: the file contained chromosome names <str>
    -o | --output: the output folder <str>
    -p | --prefix: the prefix of output <str>
    -t | --thread: # of cpus <int>
    --cngain: cutoff for CN gain considered by AA [5] <float>
    --cnsize: cutoff for CN interval size considered by AA [50000] <int>
    --paired: the input bam is paired-end <bool>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o ha:b:c:o:p:t:, --long anno:,bam:,chrom:,output:,prefix:,cngain:,cnsize:,thread:, \
  --long help,paired, \
  -- "$@"`


if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREAD=10
ANNO_DIR=
BAM_FILE=
CHROM_FILE=
PREFIX="cnvReadDepth"
CN_GAIN=5
CN_SIZE=50000
OUTPUT_DIR="./"
PAIRED_FLAG=false

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -a | --anno ) ANNO_DIR="$2"; shift 2 ;;
    -b | --bam ) BAM_FILE="$2"; shift 2 ;;
    -c | --chrom ) CHROM_FILE="$2"; shift 2 ;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2 ;;
    -p | --prefix ) PREFIX="$2"; shift 2 ;;
    -t | --thread ) THREAD="$2"; shift 2 ;;
    --cngain ) CN_GAIN="$2"; shift 2 ;;
    --cnsize ) CN_SIZE="$2"; shift 2 ;;
    --paired ) PAIRED_FLAG=true; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

## getopt end

## check required arguments
if [ -z $BAM_FILE ]; then
  echo "-b|bam: Wrong! bam not set!"
  exit 2
else
  if [[ ! -f "$BAM_FILE" ]]; then
    echo "-b|bam: Wrong! bam not set!"
    exit 2
  fi
fi

if [[ ! -d "$ANNO_DIR" ]]; then
  echo "-a|anno: Wrong! readDepth annotations not set!"
  exit 2
fi

if [[ ! -f $CHROM_FILE ]]; then
  echo "-c|--chrom: Wrong! chrom file NOT found!"
  exit 2
fi

echo "Running pipleline with following parameters:"
echo "ANNO_DIR=$ANNO_DIR"
echo "BAM_FILE=$BAM_FILE"
echo "PREFIX=$PREFIX"
echo "CHROM_FILE=$CHROM_FILE"
echo "CN_GAIN=$CN_GAIN"
echo "CN_SIZE=$CN_SIZE"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo "PAIRED_FLAG=$PAIRED_FLAG"
echo ""


BED_RES="$OUTPUT_DIR/reads"
OUTPUT_RES="$OUTPUT_DIR/output"
PARAM_FILE="$OUTPUT_DIR/params"

# delete old files
rm -rf $OUTPUT_DIR

# make directories
if [[ ! -d $BED_RES ]]; then
  mkdir -p $BED_RES
fi

if [[ ! -d $OUTPUT_RES ]]; then
  mkdir -p $OUTPUT_RES
fi

cd $OUTPUT_DIR

echo "Generating reads from BAM..."
# preparing annotations
cp -r $ANNO_DIR ./annotations

echo "
readLength 100
fdr 0.05
overDispersion 1
gcWindowSize 100
percCNGain 0.05
percCNLoss 0.05
chunkSize 1e7
maxCores 16
readCores 16
verbose TRUE
" > $PARAM_FILE

#chrom_arr=( "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10")
#chrom_arr+=( "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20")
#chrom_arr+=( "chr21" "chr22" "chrM" "chrX" "chrY" )


echo "Generating reads from BAM..."
# extract reads from BAM per chromosome
tmpfifo="/tmp/$RANDOM.fifo"
mkfifo tmpfifo
exec 9<>tmpfifo
rm -rf tmpfifo

for((n=1;n<=${THREAD};n++))
do
  echo >&9
done

while IFS=$'\t' read -r -a chrom_arr
do
  read -u9
  ## Collapse exact duplicates
  {
    chr="${chrom_arr[0]}"
    if [ "$PAIRED_FLAG" = true ]; then
      samtools view -bf 0x40 $BAM_FILE $chr | \
        bedtools bamtobed -i stdin | \
        awk 'BEGIN{OFS="\t";FS="\t"}{$2=$2+1; print $1, $2, $3}' | \
        sort -t $'\t' -k 1,1 -k2,2n > $BED_RES/$chr.bed
    else
      samtools view -b $BAM_FILE $chr | \
        bedtools bamtobed -i stdin | \
        awk 'BEGIN{OFS="\t";FS="\t"}{$2=$2+1; print $1, $2, $3}' | \
        sort -t $'\t' -k 1,1 -k2,2n > $BED_RES/$chr.bed
    fi
    echo >&9
  } &
done < $CHROM_FILE
wait
exec 9>&-
exec 9<&-

echo "Reads done!"

echo "Running ReadDepth..."
# running ReadDepth
readDepthRun.R --input ./ > readDepth.log 2>&1

echo "ReadDepth done!"

cp $OUTPUT_RES/alts.dat ${PREFIX}.CNVs.bed

echo "filtering CNVs, span_with >= ${CN_SIZE}bp, copy_number >= ${CN_GAIN}!"
# filtering CNVs, span_with >= ${CN_SIZE}kb, copy_number >= ${CN_GAIN}
awk -v gain="$CN_GAIN" -v size="$CN_SIZE" 'BEGIN{OFS="\t";FS="\t";}
{
  if ($5 >= gain && ($3-$2) >= size) {
    print
  }
}' ${PREFIX}.CNVs.bed > ${PREFIX}.CNVs_GAIN.bed

rm -rf ./annotations
rm -rf $BED_RES

echo "cnvReadDepth.sh done!"
