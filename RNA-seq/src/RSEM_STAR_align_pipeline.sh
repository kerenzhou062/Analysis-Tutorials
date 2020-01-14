#!/bin/bash
#SBATCH --job-name=RSEM_STAR_align_pipeline    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 10                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=60G                      # Amount of memory in GB
#SBATCH --time=120:10:00               # Time limit hrs:min:sec
#SBATCH --output=RSEM_STAR_align_pipeline.log               # Time limit hrs:min:sec

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
# ref:https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/DAC/STAR_RSEM.sh

# output: all in the working directory, fixed names
# ${sortedGenomeBAM}                 # alignments, standard sorted BAM, agreed upon formatting
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Quant.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
# Quant.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
# Quant.pdf                                     # RSEM diagnostic plots
# Signal.{Unique,UniqueMultiple}.strand{plus,minus}.rpm.bw # 4 bigWig files for stranded data
# Signal.{Unique,UniqueMultiple}.unstranded.rpm.bw  # 2 bigWig files for unstranded data

function showHelp {
  echo -ne "usage: sbatch RSEM_STAR_align_pipeline.sh -n <thread_num> -o <log> --mem <30G> "
  echo -ne "RSEM_STAR_align_pipeline.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -o | --output: the output folder <str>
    -p | --prefix: the output prefix (align) <str>
    -t | --thread: # of cpus <int>
    --max-mismatch: --outFilterMismatchNoverLmax in STAR (0.04) <float>
    --mem: # of memory used for RSEM --ci-memory [GB] (30) <int>
    --rsem-genome-dir: RSEM genome directory - prepared with RSEM_STAR_prep.sh <str>
    --star-genome-dir: STAR genome directory - prepared with RSEM_STAR_prep.sh <str>
    --read1: fastq of read1 <str>
    --read2: fastq of read2 (not set if single-end) <str>
    --strandedess: strandedness of the RNA-Seq reads ('none[default]|forward|reverse') <str>
    --seq-type: RNA-seq type, possible values: PE SE
    --append-names: set RSEM with --append-names <bool>
    --disable-bw: do not generate bigWig files <str>
    --skip-mapping: skip reads mapping step <bool>
    --skip-txsort: skip toTranscriptome.out.bam sorting step <bool>
    --zcat-flag: set if the input fastq <bool>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o hp:o:t: --long help,thread:,max-mismatch:,mem:,prefix: \
  --long rsem-genome-dir:,star-genome-dir:,read1:,read2:,strandedness:,seq-type:, \
  --long append-names,disable-bw,skip-mapping,skip-txsort,zcat-flag \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREAD=10
MEMORY="60"
PREFIX="align"
OUTPUT_DIR=
STAR_GENOME_DIR=
RSEM_GENOME_DIR=
MAX_MISMATCH=0.04
SEQ_TYPE="str_PE"
READ1=
READ2=""
STRANDEDNESS="none"
APPEND_NAMES=false
DISABLE_BW=false
ZCAT_FLAG=false
SKIP_MAPPING=false
SKIP_TX_SORT=false

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2 ;;
    -p | --prefix ) PREFIX="$2"; shift 2 ;;
    -t | --thread ) THREAD="$2"; shift 2 ;;
    --rsem-genome-dir ) RSEM_GENOME_DIR="$2"; shift 2 ;;
    --star-genome-dir ) STAR_GENOME_DIR="$2"; shift 2 ;;
    --max-mismatch ) MAX_MISMATCH="$2"; shift 2 ;;
    --mem ) MEMORY="$2"; shift 2 ;;
    --read1 ) READ1="$2"; shift 2 ;;
    --read2 ) READ2="$2"; shift 2 ;;
    --strandedness ) STRANDEDNESS="$2"; shift 2 ;;
    --seq-type ) SEQ_TYPE="$2"; shift 2 ;;
    --append-names ) APPEND_NAMES=true; shift ;;
    --disable-bw ) DISABLE_BW=true; shift ;;
    --skip-mapping ) SKIP_MAPPING=true; shift ;;
    --skip-txsort ) SKIP_TX_SORT=true; shift ;;
    --zcat-flag ) ZCAT_FLAG=true; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

## getopt end

## check required arguments
if [ -z $READ1 ]; then
  echo "--read1: Wrong! read1 NOT found!"
  exit 2
fi

case "$SEQ_TYPE" in
PE|SE)
  #Nothing: stranded data
  ;;
*)
  echo "--seq-type: Wrong! Possible values: PE SE!"
  exit 2
  ;;
esac

case "$STRANDEDNESS" in
none|forward|reverse)
  #Nothing: stranded data
  ;;
*)
  echo "--strandedness: Wrong! Possible values: none forward reverse!"
  exit 2
  ;;
esac

if [[ $SEQ_TYPE == "SE" ]]; then
  READ2=""
else
  if [ -z $READ2 ]; then
    echo "--read2: Wrong! read2 NOT found!"
    exit 2
  fi
fi

echo "Running pipleline with following parameters:"
echo "THREAD=$THREAD"
echo "MEMORY=$MEMORY"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo "PREFIX=$PREFIX"
echo "STAR_GENOME_DIR=$STAR_GENOME_DIR"
echo "RSEM_GENOME_DIR=$RSEM_GENOME_DIR"
echo "MAX_MISMATCH=$MAX_MISMATCH"
echo "SEQ_TYPE=$SEQ_TYPE"
echo "READ1=$READ1"
echo "READ2=$READ2"
echo "APPEND_NAMES=$APPEND_NAMES"
echo "STRANDEDNESS=$STRANDEDNESS"
echo "DISABLE_BW=$DISABLE_BW"
echo "SKIP_MAPPING=$SKIP_MAPPING"
echo "SKIP_TX_SORT=$SKIP_TX_SORT"
echo "ZCAT_FLAG=$ZCAT_FLAG"
echo ""

# executables
STAR=STAR                             
RSEM=rsem-calculate-expression        
bedGraphToBigWig=bedGraphToBigWig

zcatCommand=""
if $ZCAT_FLAG; then
  zcatCommand="--readFilesCommand zcat"
fi

if [[ ! -d $OUTPUT_DIR ]]; then
  mkdir -p $OUTPUT_DIR
else
  rm -rf $OUTPUT_DIR
  mkdir -p $OUTPUT_DIR
fi

cd $OUTPUT_DIR

# to avoid out of memory 
MEMORY=$((MEMORY-3))
# STAR parameters: common
STARparCommon=" --genomeDir $STAR_GENOME_DIR  --readFilesIn ${READ1} ${READ2}   --outSAMunmapped Within --outFilterType BySJout \
  --outSAMattributes NH HI AS NM MD    --outFilterMultimapNmax 20   --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax ${MAX_MISMATCH}   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000 \
  --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 $zcatCommand "

# STAR parameters: run-time, controlled by DCC
STARparRun=" --runThreadN ${THREAD} --genomeLoad NoSharedMemory  --limitBAMsortRAM 0"

# STAR parameters: type of BAM output: quantification or sorted BAM or both
#     OPTION: sorted BAM output
## STARparBAM="--outSAMtype BAM SortedByCoordinate"
#     OPTION: transcritomic BAM for quantification
## STARparBAM="--outSAMtype None --quantMode TranscriptomeSAM"
#     OPTION: both
STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM "

# STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 

if [[ $STRANDEDNESS == "none" ]]; then
  STARparStrand="--outSAMstrandField intronMotif"
  STARparWig="--outWigStrand Unstranded"
else
  STARparStrand=""
  STARparWig="--outWigStrand Stranded"
fi

# STAR parameters: metadata
STARparsMeta="--outSAMheaderCommentFile commentsENCODElong.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate"

## not needed ## --outSAMheaderPG @PG ID:Samtools PN:Samtools CL:"$samtoolsCommand" PP:STAR VN:0.1.18"
ANNO_NAME=$(basename $STAR_GENOME_DIR)
# ENCODE metadata BAM comments
echo -e '@CO\tLIBID:ENCLB175ZZZ
@CO\tREFID:ENCFF001RGS
@CO\tANNID:custom.annotation.${ANNO_NAME}
@CO\tSPIKEID:ENCFF001RTP VN:Ambion-ERCC Mix, Cat no. 445670' > commentsENCODElong.txt

# rename bam
sortedGenomeBAM="${PREFIX}.sortedByCoord.out.bam"
txBAM="${PREFIX}.toTranscriptome.out.bam"

if $SKIP_MAPPING; then
  echo "Skip mapping."
else
  ###### STAR command
  echo $STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta
  $STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta
  echo "Finish mapping."
  echo "Rename output bam..."
  mv Aligned.sortedByCoord.out.bam ${sortedGenomeBAM}
  mv Aligned.toTranscriptome.out.bam ${txBAM}
  echo -e "index ${sortedGenomeBAM}..."
  samtools index -b ${sortedGenomeBAM}
fi

if $DISABLE_BW; then
  echo "Skip signal tracks generation."
else
  echo "Generating bedGraph signal tracks with raw counts..."
  ###### bedGraph generation, now decoupled from STAR alignment step
  # working subdirectory for this STAR run
  mkdir Signal_RAW
  
  echo "$STAR --runMode inputAlignmentsFromBAM   --inputBAMfile ${sortedGenomeBAM} "
  echo "  --outWigType bedGraph $STARparWig --outWigNorm None --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr"
  $STAR --runMode inputAlignmentsFromBAM   --inputBAMfile ${sortedGenomeBAM} \
    --outWigType bedGraph $STARparWig --outWigNorm None --outFileNamePrefix ./Signal_RAW/ --outWigReferencesPrefix chr
  
  echo "Generating bedGraph signal tracks nomalized by RPM..."
  ###### bedGraph generation, now decoupled from STAR alignment step
  # working subdirectory for this STAR run
  mkdir Signal_RPM
  
  echo "$STAR --runMode inputAlignmentsFromBAM   --inputBAMfile ${sortedGenomeBAM} "
  echo "  --outWigType bedGraph $STARparWig --outWigNorm RPM --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr"
  $STAR --runMode inputAlignmentsFromBAM   --inputBAMfile ${sortedGenomeBAM} \
    --outWigType bedGraph $STARparWig --outWigNorm RPM --outFileNamePrefix ./Signal_RPM/ --outWigReferencesPrefix chr
  
  echo "Converting bedGraph tracks to bigWig tracks..."

  ###### bigWig conversion commands
  # exclude spikeins
  grep ^chr $STAR_GENOME_DIR/chrNameLength.txt > chrNL.txt
  
  if [[ $STRANDEDNESS == "none" ]]; then
    # unstranded data
    for imult in Unique UniqueMultiple
    do
      ## raw counts 
      grep ^chr ./Signal_RAW/Signal.$imult.str1.out.bg | LC_COLLATE=C sort -k1,1 -k2,2n > sig.tmp
      $bedGraphToBigWig sig.tmp chrNL.txt  Signal.$imult.unstranded.raw.bw
      ## RPM 
      grep ^chr ./Signal_RPM/Signal.$imult.str1.out.bg | LC_COLLATE=C sort -k1,1 -k2,2n > sig.tmp
      $bedGraphToBigWig sig.tmp chrNL.txt  Signal.$imult.unstranded.rpm.bw
    done
  else
    # stranded data
    if [[ STRANDEDNESS == "reverse" ]]; then
      str[1]="minus";
      str[2]="plus";
    else
      str[1]="plus";
      str[2]="minus";
    fi
    
    for istr in 1 2
    do
      for imult in Unique UniqueMultiple
      do
        ## raw counts 
        grep ^chr ./Signal_RAW/Signal.$imult.str${istr}.out.bg | LC_COLLATE=C sort -k1,1 -k2,2n > sig.tmp
        $bedGraphToBigWig sig.tmp  chrNL.txt Signal.$imult.${str[istr]}.raw.bw
        ## RPM 
        grep ^chr ./Signal_RPM/Signal.$imult.str${istr}.out.bg | LC_COLLATE=C sort -k1,1 -k2,2n > sig.tmp
        $bedGraphToBigWig sig.tmp  chrNL.txt Signal.$imult.${str[istr]}.rpm.bw
      done
    done
  fi

  echo "Rename bigWig tracks..."
  rename "Signal." "${PREFIX}." *.bw
  mv Signal_RAW/Log.out "${PREFIX}.Signal.raw.Log.out"
  mv Signal_RPM/Log.out "${PREFIX}.Signal.rpm.Log.out"

  echo "Deleting bedGraph tracks..."
  rm -rf ./Signal_RAW
  rm -rf ./Signal_RPM
fi
######### RSEM

if $SKIP_TX_SORT; then
  echo "Skip sort transcriptome BAM."
else
  #### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
  trBAMsortRAM="${MEMORY}G"
  
  echo "Sorting toTranscriptome bam..."
  mv ${txBAM} Tr.bam 
  
  if [[ $SEQ_TYPE == "SE" ]]; then
    # single-end data
    cat <( samtools view -H Tr.bam ) <( samtools view -@ ${THREAD} Tr.bam | sort -S ${trBAMsortRAM} -T ./ ) | \
      samtools view -@ ${THREAD} -bS - > ${txBAM}
  else
    # paired-end data, merge mates into one line before sorting, and un-merge after sorting
    cat <( samtools view -H Tr.bam ) <( samtools view -@ ${THREAD} Tr.bam | \
      awk '{printf "%s", $0 " "; getline; print}' | sort -S ${trBAMsortRAM} -T ./ | tr ' ' '\n' ) | \
      samtools view -@ $THREAD -bS - > ${txBAM}
  fi
  ## delete temp bam
  rm -f Tr.bam
fi

echo "Running RSEM: ${RSEM}..."
# RSEM parameters: common
## --estimate-rspd enabls RSEM to learn from data how the reads are distributed across a transcript. 
## The learned statistics can help us assess if any positional biases are shown in the data.
## The --calc-ci option tells RSEM to calculate credibility intervals and CQV values for each isoform / gene
## The coefficient of quartile variation (CQV), which is a robust way to measure 
## the ratio between standard deviation and mean.
## Small CQV(0.05) means that we have enough reads to produce a good estimate
if $APPEND_NAMES; then
  RSEMparCommon="--bam --append-names --estimate-rspd  --calc-ci --no-bam-output --seed 12345"
else
  RSEMparCommon="--bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345"
fi

mb_memory=$((MEMORY*1000))
# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun=" -p $THREAD --ci-memory ${MEMORY} "

# RSEM parameters: data type dependent

if [[ $SEQ_TYPE == "PE" ]]; then
  RSEMparType="--paired-end"
else
  RSEMparType=""
fi

if [[ $STRANDEDNESS == "none" ]]; then
  RSEMparType="$RSEMparType --strandedness none"
elif [[ $STRANDEDNESS == "reverse" ]]; then
  RSEMparType="$RSEMparType --strandedness reverse"
else
  RSEMparType="$RSEMparType --strandedness forward"
fi

###### RSEM command
echo "$RSEM $RSEMparCommon $RSEMparRun "
echo "  $RSEMparType Aligned.toTranscriptome.out.bam $RSEM_GENOME_DIR Quant >& Log.rsem"
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType ${txBAM} $RSEM_GENOME_DIR Quant >& Log.rsem

###### RSEM diagnostic plot creation
# Notes:
# 1. rsem-plot-model requires R (and the Rscript executable)
# 2. This command produces the file Quant.pdf, which contains multiple plots
echo "rsem-plot-model Quant Quant.pdf"
rsem-plot-model Quant Quant.pdf

## deleting temp files
echo "Deleting temp files..."
rm -rf _STARtmp
rm -f *out.bg
rm -f sig.tmp

echo "Rename outputs..."
#find . -type f -name "*.rpm.bw" | perl -pe 'print $_; s/Signal/HepG2_control_rep2/' | xargs -n2 mv
mv Quant.genes.results "${PREFIX}.genes.results"
mv Quant.isoforms.results "${PREFIX}.isoforms.results"
mv Quant.pdf "${PREFIX}.Quant.pdf"
mv Log.rsem "${PREFIX}.rsem.log"
mv SJ.out.tab "${PREFIX}.SJ.out.tab"

rename "Aligned." "${PREFIX}." *.out
rename "Log." "${PREFIX}.Log." *.out

echo "RSEM_STAR pipeline done."
