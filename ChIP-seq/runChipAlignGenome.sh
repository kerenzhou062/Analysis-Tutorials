#!/bin/bash

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
# ref:https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c

function showHelp {
  echo -ne "usage: runChipAlignGenome.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -f | --fasta: genome fasta file  <str>
    -i | --index: genome index <str>
    -o | --output: the output directory (eg. alignment/HepG2_CTCF_WT/IP/rep1) <str>
    -p | --prefix: the output prefix (eg. HepG2_CTCF_WT_IP_rep1) <str>
    -q | --quality: the mapping quality cutoff [30] <str>
    -r | --root: the root IP/input directory (eg. alignment/HepG2_CTCF_WT/IP) <str>
    -t | --thread: # of cpus [10] <int>
    --fastq1: fastqs of read1 (fq_R1_run1,fq_R1_run2,fq_R1_run3) <str>
    --fastq2: fastqs of read2 (fq_R2_run1,fq_R2_run2,fq_R2_run3) (not set if single-end) <str>
    --trim: trim length <int>
    --load: load necessary modules on slurm server <int>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o hf:i:p:q:o:r:t:, --long help,thread:,fastq1:,fastq2:,trim:, \
  --long fasta:,index:,output:,prefix:,quality:,root:, \
  --long load \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREADS=10
PREFIX=
MAPQ_THRESH=30
FASTQ_R1=""
FASTQ_R2=""
TRIM_LEN=50
GENOME_INDEX=
GENOME_FASTA=
ROOT_DIR=
OUTPUT_DIR=
LOAD_FLAG=false

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -f | --fasta ) GENOME_FASTA="$2"; shift 2 ;;
    -i | --index ) GENOME_INDEX="$2"; shift 2 ;;
    -p | --prefix ) PREFIX="$2"; shift 2 ;;
    -q | --quality ) MAPQ_THRESH="$2"; shift 2 ;;
    -r | --root ) ROOT_DIR="$2"; shift 2 ;;
    -t | --thread ) THREADS="$2"; shift 2 ;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2 ;;
    --fastq1 ) FASTQ_R1="$2"; shift 2 ;;
    --fastq2 ) FASTQ_R2="$2"; shift 2 ;;
    --trim ) TRIM_LEN="$2"; shift 2 ;;
    --load ) LOAD_FLAG=true; shift ;;
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

if [ -z "$GENOME_FASTA" ] || [ ! -f "$GENOME_FASTA" ]; then
  echo "--fasta: Wrong! read1 NOT found!"
  exit 2
fi

if [ -z "$ROOT_DIR" ]; then
  echo "--root: Wrong! root directory not SET!"
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

echo "Running pipleline with following parameters:"
echo "THREADS=$THREADS"
echo "PREFIX=$PREFIX"
echo "FASTQ_R1=$FASTQ_R1"
echo "FASTQ_R2=$FASTQ_R2"
echo "TRIM_LEN=$TRIM_LEN"
echo "GENOME_INDEX=$GENOME_INDEX"
echo "GENOME_FASTA=$GENOME_FASTA"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo ""

if [ "$LOAD_FLAG" = true ]; then
  echo -e "load modules..."
  module load picard/2.21.1
  module load bedtools/2.29.0
  module load Python/3.6.5
  module load Anaconda2/4.2.0
fi

# ===========step 0a==============
# public bariables
# ================================
echo -e "\nPreparing basic variables..."

# ===========step 0b==============
# cat multiple run fastqs if necessary
# ================================
## other arguments
OUTPUT_PREFIX="${OUTPUT_DIR}/${PREFIX}"
ALIGN_LOG="${OUTPUT_PREFIX}.align.log"
FLAGSTAT_QC="${OUTPUT_PREFIX}.flagstat.qc"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
else
  rm -rf "$OUTPUT_DIR"
  mkdir -p "$OUTPUT_DIR"
fi

cd  "$OUTPUT_DIR"

# ===========step 0==============
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

# ===========step 1a==============
# Read alignment Single-end (bowtie2 aligner)
# ================================
# output:
# ● BAM file $bam
# ● mapping stats from flagstat (SAMstats) $FLAGSTAT_QC
echo "Step 1a starting..."
echo -e "Start read alignment (bowtie2 aligner)..."

if [ -z "$FASTQ_R2" ]; then
  # mapping single-end reads
  zcat -f ${FASTQ_R1} | bowtie2 --mm -x ${GENOME_INDEX} \
    --threads ${THREADS} -U - 2> ${ALIGN_LOG} | \
    samtools view -Su /dev/stdin | \
    samtools sort -T ${OUTPUT_PREFIX} -O bam -o ${OUTPUT_PREFIX}.bam
else
  # mapping paired-end reads
  bowtie2 -X2000 --mm -x ${GENOME_INDEX} --threads ${THREADS}   -1 $FASTQ_R1 -2 $FASTQ_R2 2> ${ALIGN_LOG} | \
    samtools view -Su /dev/stdin | \
    samtools sort -T ${OUTPUT_PREFIX} -O bam -o ${OUTPUT_PREFIX}.bam
fi

samtools index ${OUTPUT_PREFIX}.bam

samtools flagstat ${OUTPUT_PREFIX}.bam > ${FLAGSTAT_QC}

echo "Step 1a done!"

# =============step 1b==============
# Post-alignment filtering #
# Remove unmapped, mate unmapped
# not primary alignment, reads failing platform
# Remove low MAPQ reads
# Only keep properly paired reads
# Obtain name sorted BAM file
# ================================
# output:
# ● Filtered deduped position sorted BAM and index file ${FINAL_BAM_FILE}
# ${FINAL_BAM_INDEX_FILE}
# ● Filtered deduped name sorted BAM file ${FINAL_NMSRT_BAM_FILE}
# ● Flagstat Metric for filtered BAM file ${FINAL_BAM_FILE_MAPSTATS}
# ● Duplication metrics from MarkDuplicates ${DUP_FILE_QC}
# ● Library complexity measures ${PBC_FILE_QC}
echo "Step 1b starting..."
echo "Start post-alignment filtering..."

OFPREFIX="${OUTPUT_PREFIX}"
RAW_BAM_FILE="${OFPREFIX}.bam"
FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"

if [ -z "$FASTQ_R2" ]; then
  samtools view -F 1804 -q ${MAPQ_THRESH} \
    -b ${RAW_BAM_FILE} -o ${FILT_BAM_FILE}
else
  TMP_FILT_BAM_PREFIX="${FILT_BAM_PREFIX}.tmp.nmsrt"
  TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"
  
  # Will produce name sorted BAM
  samtools view -F 1804 -f 2 -q ${MAPQ_THRESH} -u ${RAW_BAM_FILE} | \
    samtools sort -n - -T ${TMP_FILT_BAM_PREFIX} -O bam -o ${TMP_FILT_BAM_FILE}
  
  # Remove orphan reads (pair was removed)
  # and read pairs mapping to different chromosomes
  # Obtain position sorted BAM
  samtools fixmate -r ${TMP_FILT_BAM_FILE} ${OFPREFIX}.fixmate.tmp
  samtools view -F 1804 -f 2 -u ${OFPREFIX}.fixmate.tmp | \
    samtools sort - -T ${FILT_BAM_PREFIX} -O bam -o ${FILT_BAM_FILE}
  rm ${OFPREFIX}.fixmate.tmp
  rm ${TMP_FILT_BAM_FILE}
fi

# ======================================
# Mark duplicates
# ====================================
echo "Mark duplicates..."

TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
MARKDUP="/opt/picard/2.21.1/picard.jar MarkDuplicates"
DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc" # QC file

java -Xmx4G -jar ${MARKDUP} INPUT=${FILT_BAM_FILE} \
  OUTPUT=${TMP_FILT_BAM_FILE} \
  METRICS_FILE=${DUP_FILE_QC} \
  VALIDATION_STRINGENCY=LENIENT \
  ASSUME_SORTED=true REMOVE_DUPLICATES=false

mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}

# ==========================================
# Remove duplicates
# Index final position sorted BAM
# Create final name sorted BAM
# ==========================================
echo "Remove duplicates and index BAM..."

FINAL_BAM_PREFIX="${OFPREFIX}.filt.srt.nodup"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bam.bai" # To be stored
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

if [ -z "$FASTQ_R2" ]; then
  # Index Final BAM file
  samtools view -F 1804 -b ${FILT_BAM_FILE} -o ${FINAL_BAM_FILE}
else
  FINAL_NMSRT_BAM_PREFIX="${OFPREFIX}.filt.nmsrt.nodup"
  FINAL_NMSRT_BAM_FILE="${FINAL_NMSRT_BAM_PREFIX}.bam" # To be stored

  samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} -o ${FINAL_BAM_FILE}
  samtools sort -n ${FINAL_BAM_FILE} -T ${FINAL_NMSRT_BAM_PREFIX} -o ${FINAL_NMSRT_BAM_FILE}
fi
# Index Final BAM file
samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}
samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}

# ================================
# Compute library complexity
# ================================
# for single-end
#  sort by position and strand
#  Obtain unique count statistics
# for paired-end
#  Sort by name
#  convert to bedPE and obtain fragment coordinates
#  sort by position and strand
#  Obtain unique count statistics
echo "Compute library complexity..."

# PBC File output
PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

if [ -z "$FASTQ_R2" ]; then
  # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab]
  # PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
  bedtools bamtobed -i ${FILT_BAM_FILE} | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | \
    grep -v 'chrM' | sort | uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1}
      END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
    > ${PBC_FILE_QC}
else
  # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab]
  # NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
  samtools sort -n ${FILT_BAM_FILE} -o ${OFPREFIX}.srt.tmp.bam
  bedtools bamtobed -bedpe -i ${OFPREFIX}.srt.tmp.bam | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}
      ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} 
      END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
    > ${PBC_FILE_QC}
  rm ${OFPREFIX}.srt.tmp.bam
fi
# do not delete ${FILT_BAM_FILE}
## rm ${FILT_BAM_FILE}

# link ${FINAL_BAM_FILE} to ${ROOT_DIR}
ln -f ${FINAL_BAM_FILE} ${ROOT_DIR}

echo "Step 1b done!"

# =======step 2a=================
# Convert SE & PE BAM to tagAlign (BED 3+3 format) #
# Create tagAlign file
# Create BEDPE file (for PE dataset)
# =================================
# output:
# ● tagAlign file ${FINAL_TA_FILE}
echo "Create tagAlign (& BEDPE) file..."
echo "Step 2a starting..."

FINAL_TA_FILE="${FINAL_BAM_PREFIX}.final.tagAlign.gz"
if [ -z "$FASTQ_R2" ]; then
  bedtools bamtobed -i ${FINAL_BAM_FILE} |   awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | \
    gzip -nc > ${FINAL_TA_FILE}
else
  FINAL_BEDPE_FILE="${FINAL_NMSRT_BAM_PREFIX}.bedpe.gz"

  bedtools bamtobed -bedpe -mate1 -i ${FINAL_NMSRT_BAM_FILE} | \
    awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | \
    gzip -nc > ${FINAL_TA_FILE}
fi

# link ${FINAL_TA_FILE} to ${ROOT_DIR}
ln -f ${FINAL_TA_FILE} ${ROOT_DIR}


echo "Step 2a done!"

# ===========step 2b==============
#Calculate Cross-correlation QC scores
# =======================================
# output:
# ● outFile containing NSC/RSC results in tab-delimited file of 11 columns (same file can
# be appended to from multiple runs) ${CC_SCORES_FILE}
# ● cross-correlation plot ${CC_PLOT_FILE}
echo "Step 2b starting..."
echo "Calculating Cross-correlation QC scores..."

#### for both PE and SE samples
# Trim R1 fastq to 50bp
echo "Starting trimming fastq..."
TRIM_TEMP_FASTQ="${OUTPUT_PREFIX}.trim.tmp.fastq"
TRIM_OFPREFIX="${OUTPUT_PREFIX}.trim"
TRIMMED_FASTQ_R1="${TRIM_OFPREFIX}.fastq.gz"
zcat -f ${FASTQ_R1} > ${TRIM_TEMP_FASTQ}
python $(which trimfastq.py) ${TRIM_TEMP_FASTQ} ${TRIM_LEN} | gzip -nc > ${TRIMMED_FASTQ_R1}

# Align $TRIMMED_FASTQ_R1 (not paired) with bowtie2 (step 1a SE) and use it for filtering
# step (1b) and then get $FILT_BAM_FILE (not the deduped $FINAL_BAM_FILE), which is
# filtered but not deduped.
echo "Start trimmed-reads mapping..."

TRIMLOG="${TRIM_OFPREFIX}.align.log"
TRIM_RAW_BAM_FILE="${TRIM_OFPREFIX}.bam"
TRIM_FILT_BAM_PREFIX="${TRIM_OFPREFIX}.filt.srt"
TRIM_FILT_BAM_FILE="${TRIM_FILT_BAM_PREFIX}.bam"

zcat -f ${TRIMMED_FASTQ_R1} | bowtie2 --mm -x ${GENOME_INDEX} --threads ${THREADS} -U - 2> ${TRIMLOG} | \
  samtools view -Su /dev/stdin | \
  samtools sort -T ${TRIM_OFPREFIX} -O bam -o ${TRIM_RAW_BAM_FILE}
samtools view -F 1804 -q ${MAPQ_THRESH} -b ${TRIM_RAW_BAM_FILE} -o ${TRIM_FILT_BAM_FILE}

# =================================
# make tagAlign for filtered (but not deduped) BAM
# and subsample it for cross-correlation analysis
# ================================
echo "Make tagAlign for trim-filtered (but not deduped) BAM..."

TEMP_TRIM_TA_FILE="${TRIM_FILT_BAM_PREFIX}.tagAlign.gz"
SUBSAMPLED_TA_FILE="${TRIM_OFPREFIX}.filt.sample.15.tagAlign.gz"

bedtools bamtobed -i ${TRIM_FILT_BAM_FILE} | \
  awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | \
  gzip -nc > ${TEMP_TRIM_TA_FILE}
zcat -f ${TEMP_TRIM_TA_FILE} | grep -v "chrM" |   shuf -n ${NREADS} --random-source=<(openssl enc -aes-256-ctr -pass \
  pass:$(zcat -f ${TEMP_TRIM_TA_FILE} | wc -c) -nosalt </dev/zero 2>/dev/null) | \
  gzip -nc > ${SUBSAMPLED_TA_FILE}

### cross-correlation analysis
# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak
# <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab>
# relPhantomPeakCoef <tab> QualityTag

### temporary block
CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"
NTHREADS=10

Rscript $(which run_spp.R) -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM \
  -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE}
sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
mv temp ${CC_SCORES_FILE}

## delete trim-related files
echo "Delete trim-related files (fastq, alignment)"
rm -f ${TRIM_TEMP_FASTQ}
rm -f ${TRIMMED_FASTQ_R1}
rm -f ${TRIM_RAW_BAM_FILE}
rm -f ${TRIM_FILT_BAM_FILE}
rm -f ${TEMP_TRIM_TA_FILE}
rm -f ${SUBSAMPLED_TA_FILE}

echo "Step 2b done!"

# ==========step 2c==================
# Create pseudoReplicates
# input TagAlign file $FINAL_TA_FILE
# =====================================
# output:
# 2 pseudoreplicate virtual SE tagAlign files ${PR1_TA_FILE} ${PR2_TA_FILE}
echo "Step 2c starting..."
echo "Creating pseudoReplicates..."


PR_PREFIX="${OFPREFIX}.filt.nodup"
PR1_TA_FILE="${PR_PREFIX}.pr1.tagAlign.gz"
PR2_TA_FILE="${PR_PREFIX}.pr2.tagAlign.gz"

if [ -z "$FASTQ_R2" ]; then
  # Get total number of read pairs
  nlines=$( zcat -f ${FINAL_TA_FILE} | wc -l )
  nlines=$(( (nlines + 1) / 2 ))
  # Shuffle and split BED file into 2 equal parts
  zcat -f ${FINAL_TA_FILE} | shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f \
    ${FINAL_TA_FILE} | wc -c) -nosalt </dev/zero 2>/dev/null) | split -d -l ${nlines} - ${PR_PREFIX}
  # Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01
  # Convert reads into standard tagAlign file
  gzip -nc "${PR_PREFIX}00" > ${PR1_TA_FILE}
  rm "${PR_PREFIX}00"
  gzip -nc "${PR_PREFIX}01" > ${PR2_TA_FILE}
  rm "${PR_PREFIX}01"
else
  joined="${FINAL_BAM_PREFIX}.final.tagAlign.temp.bedpe"
  # Make temporary fake BEDPE file from FINAL_TA_FILE
  zcat -f ${FINAL_TA_FILE} | sed 'N;s/\n/\t/' > $joined
  # Get total number of read pairs
  nlines=$( zcat -f ${joined} | wc -l )
  nlines=$(( (nlines + 1) / 2 ))
  # Shuffle and split BEDPE file into 2 equal parts
  zcat -f ${joined} | shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f \
    ${FINAL_TA_FILE} | wc -c) -nosalt </dev/zero 2>/dev/null) | split -d -l ${nlines} - ${PR_PREFIX}
  # Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01
  # Convert fake BEDPE into standard tagAlign file
  awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",
    $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "${PR_PREFIX}00" | gzip -nc > ${PR1_TA_FILE}
  rm "${PR_PREFIX}00"
  awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",
    $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "${PR_PREFIX}01" | gzip -nc > ${PR2_TA_FILE}
  rm "${PR_PREFIX}01"
  rm -f ${joined}
fi

# link ${PR1_TA_FILE} and ${PR2_TA_FILE} to ${ROOT_DIR}
ln -f ${PR1_TA_FILE} ${ROOT_DIR}
ln -f ${PR2_TA_FILE} ${ROOT_DIR}

echo "Step 2c done!"

echo "Skip step to 2g..."
# ========step 2g=============
# Calculate GC bias
# =======================
# output:
# ● GC bias Plot ${GC_BIAS_PLOT}
# ● GC bias log ${GC_BIAS_LOG}
echo "Step 2g starting..."

# we don't use plot directly generated from picard
# we process picard’s text output and make a plot
# FINAL_BAM_PREFIX="${OFPREFIX}.filt.srt.nodup"
# FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored

GCBIAS="/opt/picard/2.21.1/picard.jar CollectGcBiasMetrics"
GC_BIAS_LOG="${FINAL_BAM_PREFIX}.gc_bias.txt"
GC_BIAS_PLOT="${FINAL_BAM_PREFIX}.gc_bias.pdf"
SUMMARY="${FINAL_BAM_PREFIX}.gc_bias.summary.txt"

java -Xmx6G -XX:ParallelGCThreads=1 -jar \
  ${GCBIAS} R="${GENOME_FASTA}" I="${FINAL_BAM_FILE}" O="${GC_BIAS_LOG}" \
  USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE \
  VERBOSITY=ERROR QUIET=TRUE \
  ASSUME_SORTED=FALSE \
  CHART="${GC_BIAS_PLOT}" S="${SUMMARY}"
# use ${GC_BIAS_LOG} into the following pyhton script
# data_file: ${GC_BIAS_LOG}
# prefix: any good prefix for output file name

source activate py3
plotGC.py -data "${GC_BIAS_LOG}" -prefix "${FINAL_BAM_PREFIX}"

echo "Step 2g done!"

if [ "$multipleR1" = true ]; then
  echo "Delete $FASTQ_R1"
  rm -rf "$FASTQ_R1"
fi

if [ "$multipleR2" = true ]; then
  echo "Delete $FASTQ_R2"
  rm -rf "$FASTQ_R2"
fi
