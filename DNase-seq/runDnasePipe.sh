#!/bin/bash

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
# ref:https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c

function showHelp {
  echo -ne "usage: runDnasePipe.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -a | --adapter: the adapter file ({barcode}_5P<tab>seq), work with paired-end data <str>
    -b | --barcode: the barcode name prefix in adapter file ({barcode}_5P and {barcode}_7P), work with paired-end data <str>
    -c | --chrsize: the chromsome size file <str>
    -g | --genome: the genome [hg19|hg38] (hg19) <str>
    -i | --index: the BWA genome index <str>
    -m | --mappable: the tgz file of mappable files (mappable_target.starch, center_sites.starch, chrom_sizes.bed) <str>
    -n | --samsize: the sample size (reads number) for running qc [15000000] <int>
    -o | --output: the output directory <str>
    -p | --prefix: the output prefix <str>
    -q | --quality: the mapping quality cutoff [10] <int>
    -r | --rlen: the read length [36] <int>
    -s | --seqtype: the sequencing library (PE|SE) <str>
    -t | --thread: # of cpus [10] <int>
    --fastq1: fastqs of read1 (fq_R1_run1,fq_R1_run2,fq_R1_run3) <str>
    --fastq2: fastqs of read2 (fq_R2_run1,fq_R2_run2,fq_R2_run3) (not set if single-end) <str>
    --minor: also call hotspots on chrM and scaffolds <bool>
    --umi: whether reads in bam contain UMI ids <bool>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o ha:b:c:g:i:m:n:o:p:q:r:s:t, --long help,adapter:,barcode:,chrsize:,genome:, \
  --long index:,mappable:,samsize:,output:,prefix:,quality:,rlen:,seqtype:,thread:,fastq1:,fastq2:, \
  --long minor,umi \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREADS=10
ADAPTER_FILE=
BARCODE=
CHROM_SIZE_FILE=
GENOME="hg19"
GENOME_INDEX=
MAPPABLE_TGZ=
PREFIX=
MAPQ_THRESH=10
SEQ_TYPE="SE"
FASTQ_R1=""
FASTQ_R2=""
READ_LEN=36
SAMPLE_SIZE=15000000
OUTPUT_DIR=
MINOR_FLAG=false
UMI_FLAG=false

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -a | --adapter ) ADAPTER_FILE="$2"; shift 2 ;;
    -b | --barcode ) BARCODE="$2"; shift 2 ;;
    -c | --chrsize ) CHROM_SIZE_FILE="$2"; shift 2 ;;
    -g | --genome ) GENOME="$2"; shift 2 ;;
    -i | --index ) GENOME_INDEX="$2"; shift 2 ;;
    -m | --mappable ) MAPPABLE_TGZ="$2"; shift 2 ;;
    -n | --samsize ) SAMPLE_SIZE="$2"; shift 2 ;;
    -p | --prefix ) PREFIX="$2"; shift 2 ;;
    -q | --quality ) MAPQ_THRESH="$2"; shift 2 ;;
    -r | --rlen ) READ_LEN="$2"; shift 2 ;;
    -s | --seqtype ) SEQ_TYPE="$2"; shift 2 ;;
    -t | --thread ) THREADS="$2"; shift 2 ;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2 ;;
    --fastq1 ) FASTQ_R1="$2"; shift 2 ;;
    --fastq2 ) FASTQ_R2="$2"; shift 2 ;;
    --minor ) MINOR_FLAG=true; shift ;;
    --umi ) UMI_FLAG=true; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done
## getopt end

## check required arguments
if [ -z "$FASTQ_R1" ]; then
  echo "--fastq1: Wrong! read1 NOT found!"
  exit 2
fi

if [ -z "$CHROM_SIZE_FILE" ] || [ ! -f "$CHROM_SIZE_FILE" ]; then
  echo "--chrsize: Wrong! chromsome size file NOT found!"
  exit 2
fi

if [ -z "$MAPPABLE_TGZ" ]; then
  echo "--mappable: Wrong! mappable tgz file not SET!"
  exit 2
fi

if [ -z "$OUTPUT_DIR" ]; then
  echo "--output: Wrong! output directory not SET!"
  exit 2
fi

if [ -z "$PREFIX" ]; then
  echo "--prefix: Wrong! prefix for output not SET!"
  exit 2
fi

if [ "$GENOME" != "hg19" ] && [ "$GENOME" != "h38" ]; then
  echo "--genome: Wrong! Should be hg19 or hg38!"
  exit 2
fi

if [ "$SEQ_TYPE" != "SE" ] && [ "$SEQ_TYPE" != "PE" ]; then
  echo "--seqtype: Wrong! Should be SE or PE!"
  exit 2
fi

echo "Running pipleline with following parameters:"
echo "THREADS=$THREADS"
echo "ADAPTER_FILE=$ADAPTER_FILE"
echo "BARCODE=$BARCODE"
echo "CHROM_SIZE_FILE=$CHROM_SIZE_FILE"
echo "GENOME=$GENOME"
echo "GENOME_INDEX=$GENOME_INDEX"
echo "MAPPABLE_TGZ=$MAPPABLE_TGZ"
echo "PREFIX=$PREFIX"
echo "MAPQ_THRESH=$MAPQ_THRESH"
echo "SEQ_TYPE=$SEQ_TYPE"
echo "FASTQ_R1=$FASTQ_R1"
echo "FASTQ_R2=$FASTQ_R2"
echo "READ_LEN=$READ_LEN"
echo "SAMPLE_SIZE=$SAMPLE_SIZE"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo "MINOR_FLAG=$MINOR_FLAG"
echo "UMI_FLAG=$UMI_FLAG"
echo ""

# ===========step 0a==============
# public bariables
# copy required files to output
# ================================
echo -e "\nPreparing basic variables..."

if [[ "$SEQ_TYPE" == "SE" ]]; then
  SEQ_TYPE="se"
else
  SEQ_TYPE="pe"
fi

# change direcotry to output

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
else
  rm -rf "$OUTPUT_DIR"
  mkdir -p "$OUTPUT_DIR"
fi

cd $OUTPUT_DIR

# copy required files to output

cp $MAPPABLE_TGZ $OUTPUT_DIR

MAPPABLE_TGZ=$(basename $MAPPABLE_TGZ)
MAPPABLE_TGZ="$OUTPUT_DIR/$MAPPABLE_TGZ"

if [[ ! -z "$ADAPTER_FILE" ]]; then
  cp $ADAPTER_FILE $OUTPUT_DIR
  ADAPTER_FILE=$(basename $ADAPTER_FILE)
  ADAPTER_FILE="$OUTPUT_DIR/$ADAPTER_FILE"
fi

# ===========step 1==============
# Try to cat multiple runs for 
# input fastqs
# ================================
multipleR1=false
multipleR2=false

IFS_OLD=$IFS
IFS=',' read -d '' -ra fqArr1 < <(printf '%s\0' "$FASTQ_R1")
IFS=$IFS_OLD

fqArr1Length="${#fqArr1[@]}"

if [ "$fqArr1Length" -gt "1" ]; then
  multipleR1=true
  FASTQ_R1="${OUTPUT_DIR}/${OUTPUT_PREFIX}_R1.fastq"
  catCommand="zcat -f"
  echo "Cating fastq files for multiple runs of --fastq1..."
  for i in "${!fqArr1[@]}"
  do
    fq=$(realpath ${fqArr1[$i]})
    if [ ! -f $fq ]; then
      echo "--fastq1: Wrong! Reads '$fq' NOT found!"
      exit 2
    else
      catCommand="${catCommand} $fq"
    fi
  done
  catCommand="$catCommand > ${FASTQ_R1}"
  echo "$catCommand"
  $catCommand
else
  if [ ! -f "$FASTQ_R1" ]; then
    echo "--fastq1: Wrong! read1 NOT found!"
    exit 2
  fi
fi

if [ ! -z "$FASTQ_R2" ]; then
  IFS_OLD=$IFS
  IFS=',' read -d '' -ra fqArr2 < <(printf '%s\0' "$FASTQ_R2")
  IFS=$IFS_OLD
  fqArr2Length="${#fqArr2[@]}"
  if [ "$fqArr2Length" -gt "1" ]; then
    multipleR2=true
    FASTQ_R2="${OUTPUT_DIR}/${OUTPUT_PREFIX}_R2.fastq"
    catCommand="zcat -f"
    echo "Cating fastq files for multiple runs of --fastq2..."
    for i in "${!fqArr2[@]}"
    do
      fq=$(realpath ${fqArr2[$i]})
      if [ ! -f $fq ]; then
        echo "--fastq2: Wrong! Reads '$fq' NOT found!"
        exit 2
      else
        catCommand="${catCommand} $fq"
      fi
    done
    catCommand="$catCommand > ${FASTQ_R2}"
    echo "$catCommand"
    $catCommand
  else
    if [ ! -f "$FASTQ_R2" ]; then
      echo "--fastq2: Wrong! read2 NOT found!"
      exit 2
    fi
  fi
fi

## other arguments
OUTPUT_PREFIX="${OUTPUT_DIR}/${PREFIX}"
ALIGN_LOG="${OUTPUT_PREFIX}.align.log"
FLAGSTAT_QC="${OUTPUT_PREFIX}.flagstat.qc"

pipepath=$(dirname $(realpath $(which runDnasePipe.sh)))
echo -e "Running dnase-align-bwa..."

BAM_PREFIX="$OUTPUT_DIR/$PREFIX"
FILT_BAM_PREFIX="$OUTPUT_DIR/$PREFIX.filt"
SAMPLE_SIZE=15000000

export PATH="$PATH:$pipepath/encode/dnase-align-bwa"
export PATH="$PATH:$pipepath/encode/dnase-filter"
export PATH="$PATH:$pipepath/encode/dnase-qc-bam"
export PATH="$PATH:$pipepath/encode/dnase-call-hotspots"

if [ "$SEQ_TYPE" == "se" ]; then
  # mapping single-end reads
  echo "Align reads to the genome..."
  dnase_align_bwa="$pipepath/encode/dnase-align-bwa/dnase_align_bwa_se.sh"
  $dnase_align_bwa $GENOME_INDEX $FASTQ_R1 $THREADS $BAM_PREFIX
  # filtering bam
  echo "Filter the bam..."
  dnase_filter="$pipepath/encode/dnase-filter/dnase_filter_se.sh"
  $dnase_filter ${BAM_PREFIX}.bam $MAPQ_THRESH $THREADS $FILT_BAM_PREFIX
else
  # mapping paired-end reads
  echo "Align reads to the genome..."
  dnase_align_bwa="$pipepath/dnase-align-bwa/dnase_align_bwa_pe.sh"
  $dnase_align_bwa $GENOME_INDEX $FASTQ_R1 $FASTQ_R2 $BARCODE $UMI_FLAG $ADAPTER_FILE $THREADS $BAM_PREFIX
  # filtering bam
  echo "Filter the bam..."
  dnase_filter="$pipepath/encode/dnase-filter/dnase_filter_pe.sh"
  $dnase_filter "$BAM_PREFIX.bam" $MAPQ_THRESH $THREADS $UMI_FLAG $FILT_BAM_PREFIX
fi

# running qc
HOTSPOT_DIR="$OUTPUT_DIR/"
if [[ ! -d "$HOTSPOT_DIR/hotspot-distr/data" ]]; then
  mkdir -p "$HOTSPOT_DIR/hotspot-distr/data"
fi

#echo "Quality check..."
dnase_qc_bam="$pipepath/encode/dnase-qc-bam/dnase_qc_bam.sh"
$dnase_qc_bam "$FILT_BAM_PREFIX.bam" $SAMPLE_SIZE $THREADS $SEQ_TYPE $GENOME $MAPPABLE_TGZ $READ_LEN $HOTSPOT_DIR

# running hotspot
echo "Calling broadPeak and narrowPeak..."
HOTSPOT_PREFIX="$OUTPUT_DIR/$PREFIX.broadPeak"
PEAK_PREFIX="$OUTPUT_DIR/$PREFIX.narrowPeak"
DENSITY_PREFIX="$OUTPUT_DIR/$PREFIX"
dnase_hotspot="$pipepath/encode/dnase-call-hotspots/dnase_hotspot.sh"
$dnase_hotspot "$FILT_BAM_PREFIX.bam" $CHROM_SIZE_FILE $MAPPABLE_TGZ "$HOTSPOT_PREFIX" "$PEAK_PREFIX" "$DENSITY_PREFIX" $MINOR_FLAG

# rename broadPeak and narrowPeak
mv ${HOTSPOT_PREFIX}.bed $HOTSPOT_PREFIX
mv ${PEAK_PREFIX}.bed $PEAK_PREFIX

# delete temperary files
rm -f ${BAM_PREFIX}.bam
rm -f ${BAM_PREFIX}_marked.bam
rm -f sorted.bam tmp.bam tmp_se.sam tmp.sai
rm -f ${BAM_PREFIX}.starch 
rm -f mappable_target.starch *.mappable_only.starch center_sites.starch
rm -f *.tgz

# delete cat fastqs
if [ "$multipleR1" = true ]; then
  echo "Delete $FASTQ_R1"
  rm -rf "$FASTQ_R1"
fi

if [ "$multipleR2" = true ]; then
  echo "Delete $FASTQ_R2"
  rm -rf "$FASTQ_R2"
fi

echo "Pipleline done!"
