#!/bin/bash
#SBATCH --job-name=pureCLIP_iCLIP_pipeline    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail  
#SBATCH -n 10                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=100G                      # Amount of memory in GB
#SBATCH --time=72:10:00               # Time limit hrs:min:sec
#SBATCH --output=pureCLIP_miCLIP_pipeline.log   # Standard output and error log


# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.

function showHelp {
  echo -ne "usage: sbatch pureCLIP_miCLIP_pipeline.sh -n <thread_num> -o <log> --mem <200G> "
  echo -ne "pureCLIP_miCLIP_pipeline.sh <options>\n"
  echo -e "options:
    -h | --help: show help infomation <bool>
    -b | --barcode-length: barcode length <int>
    -e | --exp-prefix: experiment prefix string <str>
    -f | --fasta: genome fasta file <str>
    -g | --gsize: genome size file <str>
    -i | --input: input fastq directory (cutadapt) <str>
    -o | --output: output result directory <str>
    -p | --pool-prefix: pooled experiment prefix <str>
    -t | --thread: # of cpus <int>
    --gtf: genome annotation gtf <str>
    --full-bed: all transcripts annotation in bed12 <str>
    --index: genome index <str>
    --longest-bed: mRNA longest annotation in bed12 <str>
    --repeat-bed: repeat bed used for filtering (eg.t/rRNA) <str>
    --skip-mapping: skip reads mapping step <bool>
    --skip-calling: skip m6A sites calling step <bool>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o hb:e:g:i:o:p:t: --long help,skip-mapping,skip-calling \
  --long input:,thread:,output:,exp-prefix:,pool-prefix:,index: \
  --long longest-bed:,full-bed:,gtf:,repeat-bed:,barcode-length:,fasta: \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREAD=10
BARCODE_LEN=9
OUTPUT_DIR=
INPUT_DIR=
EXP_PREFIX=
POOL_PREFIX=
GENOME_INDEX=
FASTA=
GTF=
GENOME_SIZE=
LONGEST_BED=
FULL_BED=
REPEAT_BED=
SKIP_MAPPING=false
SKIP_CALLING=false
while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -b | --barcode-length ) BARCODE_LEN="$2"; shift 2 ;;
    -e | --exp-prefix ) EXP_PREFIX="$2"; shift 2 ;;
    -f | --fasta ) FASTA="$2"; shift 2 ;;
    -g | --gsize ) GENOME_SIZE="$2"; shift 2 ;;
    -i | --input ) INPUT_DIR="$2"; shift 2 ;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2 ;;
    -p | --pool-prefix ) POOL_PREFIX="$2"; shift 2 ;;
    -t | --thread ) THREAD="$2"; shift 2 ;;
    --gtf ) GTF="$2"; shift 2 ;;
    --full-bed ) FULL_BED="$2"; shift 2 ;;
    --index ) GENOME_INDEX="$2"; shift 2 ;;
    --longest-bed ) LONGEST_BED="$2"; shift 2 ;;
    --repeat-bed ) REPEAT_BED="$2"; shift 2 ;;
    --skip-mapping ) SKIP_MAPPING=true; shift ;;
    --skip-calling ) SKIP_CALLING=true; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

## getopt end

## check required arguments
if [ -z $INPUT_DIR ]; then
  echo "-i|--input: Wrong! Cutadapt directory NOT found!"
  exit 2
fi

if [ -z $OUTPUT_DIR ]; then
  echo "-o|--output: Wrong! Ouput directory NOT found!"
  exit 2
fi

if [ -z $EXP_PREFIX ]; then
  echo "-e|--exp-prefix: Wrong! Experiment prefix NOT found!"
  exit 2
fi

if [ -z $POOL_PREFIX ]; then
  echo "-p|--pool-prefix: Wrong! Pooled experiment prefix NOT found!"
  exit 2
fi

if [ -z $GENOME_INDEX ]; then
  if [ -z $FASTA ]; then
    echo "-f|--fasta: Wrong! Genome fasta NOT found!"
    exit 2
  fi
  if [ -z $GTF ]; then
    echo "--gtf: Wrong! Annotation gtf NOT found!"
    exit 2
  fi
fi
# basic variables

MAP_DIR="$OUTPUT_DIR/mapping"
FASTQC_DIR="$OUTPUT_DIR/mapping/fastQC"
CROSSLINK_DIR="$OUTPUT_DIR/crosslink"
STAR_INDEX_DIR="${pureCLIP}/STAR_index"
# Read preprocessing
## Read quality $FILT_DIR
if [[ ! -d $OUTPUT_DIR ]]; then
  mkdir -p $OUTPUT_DIR
fi
## try to build index
if ! $SKIP_MAPPING; then
  #statements
  if $STAR_FLAG; then
    echo "Mapping with STAR."
    if [ -z $GENOME_INDEX ]; then
      GENOME_INDEX="$OUTPUT_DIR/STAR_index"
      if [[ ! -d $GENOME_INDEX ]]; then
        mkdir -p $GENOME_INDEX
      fi
      if [[ ! -f $GENOME_INDEX/Genome ]]; then
        ## estimate reads length
        echo "NO STAR index detected!"
        echo "Estimating reads length..."
        READ_AVG_LENGTH=100
        for i in "$INPUT_DIR/${EXP_PREFIX}*.fastq"; do
          READ_AVG_LENGTH=`awk '{if(NR%4==2) {count++; lengthSum += length} } END{print int(lengthSum/count)}' $i`
          break 1
        done
        echo "Reads length:${READ_AVG_LENGTH}"
        cd $GENOME_INDEX
        ### build index
        SJDB_OVERHANG=$((READ_AVG_LENGTH-1))
        echo "--sjdbOverhang:$SJDB_OVERHANG"
        echo "Building STAR index..."
        STAR --runThreadN ${THREAD} --runMode genomeGenerate --genomeDir $GENOME_INDEX \
        --genomeFastaFiles $FASTA --sjdbGTFfile $GTF --sjdbOverhang ${SJDB_OVERHANG} \
        -outTmpDir $MAP_DIR/_STARtmp
        echo "STAR index done!"
      else
        echo "STAR index found! Skip building index."
      fi
    else
      echo "STAR index found! Skip building index."
    fi
  fi
fi

if [[ ! -d $MAP_DIR  ]]; then
    mkdir -p $MAP_DIR
fi

if $SKIP_MAPPING ; then
  echo "Skip mapping."
else
  echo "Mapping start..."
  ## parallel start
  REP_NUM=`find $INPUT_DIR -type f -name "${EXP_PREFIX}*.trim.fastq" | wc -l`
  MAP_THREAD=$((THREAD / REP_NUM))
  
  rm -rf tmpfifo
  mkfifo tmpfifo
  exec 9<>tmpfifo
  
  for((n=1;n<=${THREAD};n++))
  do
    echo >&9
  done
  
  for i in $INPUT_DIR/${EXP_PREFIX}*.trim.fastq;
  do
    read -u9
    ## Collapse exact duplicates
    {
      cd $MAP_DIR
      PREFIX=${i%%.fastq}
      PREFIX=${PREFIX##*/} # *.trim
      MAP_PREFIX=${i%%.trim.fastq}
      MAP_PREFIX=${MAP_PREFIX##*/}
      ## collapse fastq
      fastq2collapse.pl $i - | gzip -c > ${PREFIX}.c.fastq.gz
      ## Strip random barcode (UMI)
      stripBarcode.pl -format fastq -len ${BARCODE_LEN} ${PREFIX}.c.fastq.gz - | \
        gzip -c > ${PREFIX}.c.bc.fastq.gz
      ## mapping, reduce mismatch rates
      STAR --outSAMtype BAM SortedByCoordinate --runThreadN ${MAP_THREAD} \
        --genomeDir ${GENOME_INDEX} --readFilesIn ${PREFIX}.c.bc.fastq.gz \
        --readFilesCommand  zcat --outFilterType BySJout \
        --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 --scoreDelOpen -1 \
        --alignIntronMin 20 --alignIntronMax 1000000 \
        --outFileNamePrefix ${MAP_PREFIX}. --outTmpDir ${MAP_PREFIX}_STARtmp \
        --alignEndsType EndToEnd
      rm -rf ${MAP_PREFIX}_STARtmp
      # Filtering
      ## We filter the aligned reads to obtain only reads
      ## mapping against the main chromosomes:
      samtools index -@ ${MAP_THREAD} ${MAP_PREFIX}.Aligned.sortedByCoord.out.bam
      samtools view -@ ${MAP_THREAD} -hb ${MAP_PREFIX}.Aligned.sortedByCoord.out.bam \
        -o ${MAP_PREFIX}.aligned.f.bam chr1:1 chr2:1 chr3:1 chr4:1 \
        chr5:1 chr6:1 chr7:1 chr8:1 chr9:1 chr10:1 chr11:1 chr12:1 \
        chr13:1 chr14:1 chr15:1 chr16:1 chr17:1 chr18:1 chr19:1 \
        chr20:1 chr21:1 chr22:1 chrX:1 chrY:1
      samtools view -@ ${MAP_THREAD} -hb ${MAP_PREFIX}.Aligned.sortedByCoord.out.bam \
        -o ${MAP_PREFIX}.aligned.f.bam
      samtools index -@ ${MAP_THREAD} ${MAP_PREFIX}.aligned.f.bam
      ### PCR duplicate removal using UMI
      #umi_tools dedup -I ${MAP_PREFIX}.aligned.f.bam \
      #  -S ${MAP_PREFIX}.aligned.f.duplRm.bam
      echo >&9
    } &
  done
  
  wait
  
  exec 9>&-
  exec 9<&-
  rm -rf tmpfifo
  # parallel end
fi
## Merging biological replicates

# pooling reads
cd $MAP_DIR
FIND_BAM="*.aligned.f.bam"
#MERGE_PREFIX="${POOL_PREFIX}.aligned.f.duplRm.pooled"
MERGE_PREFIX="${POOL_PREFIX}.aligned.f.pooled"
DUPLRM_BAMS=""
for i in `find ./ -maxdepth 1 -type f -name "$FIND_BAM" | sort`;
do
  DUPLRM_BAMS=`echo "${DUPLRM_BAMS}"" ""${i}"`
done

#samtools merge -@ ${MAP_THREAD} -f ${MERGE_PREFIX}.bam ${DUPLRM_BAMS}
#samtools index -@ ${MAP_THREAD} ${MERGE_PREFIX}.bam ${MERGE_PREFIX}.bam.bai

## Quality control
if [[ ! -d $FASTQC_DIR ]]; then
    mkdir $FASTQC_DIR
fi
#fastqc -o $FASTQC_DIR ${MERGE_PREFIX}.bam \
#  > ${MERGE_PREFIX}.fastqc.log 2>&1

## calling crosslink sites
if [[ ! -d $CROSSLINK_DIR  ]]; then
    mkdir -p $CROSSLINK_DIR
fi
cd $CROSSLINK_DIR

ln -sf ${MAP_DIR}/${MERGE_PREFIX}.bam.bai ./
ln -sf ${MAP_DIR}/${MERGE_PREFIX}.bam ./

# pureCLIP crosslink sites calling
if $SKIP_CALLING; then
  echo "Skip peak&site calling."
else
  echo "pureCLIP crosslink sites calling..."
  pureclip -i ${MERGE_PREFIX}.bam \
    -bai ${MERGE_PREFIX}.bam.bai \
    -g ${FASTA} -iv 'chr1;chr2;chr3;' -nt ${THREAD} \
    -o ${POOL_PREFIX}.PureCLIP.crosslink.bed \
    -or ${POOL_PREFIX}.PureCLIP.bind.bed \
    -p ${POOL_PREFIX}.PureCLIP.parameter.txt \
    > ${MERGE_PREFIX}.pureCLIP.log 2>&1
  
  echo "pureCLIP crosslink sites done."
fi
## hidden state(4.state):
## 1->0->non-enriched + non-crosslink
## 2->1->non-enriched + crosslink
## 3->2->enriched + non-crosslink
## 4->3->enriched + crosslink (green, get)
# get mutation information using CTK

#samtools fillmd -@ ${THREAD} ${MERGE_PREFIX}.bam \
#  ${FASTA} | gzip -c > ${MERGE_PREFIX}.sam.gz
#
#parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file \
#  ${POOL_PREFIX}.pooled.mutation.txt ${MERGE_PREFIX}.sam.gz ${POOL_PREFIX}.tag.bed

bedtools slop -i ${POOL_PREFIX}.PureCLIP.crosslink.bed \
  -b 5 -g ${GENOME_SIZE} | bedtools getfasta -fi ${FASTA} \
  -bed stdin -s -name+ > ${POOL_PREFIX}.PureCLIP.postive.fa

#bedtools slop -i ${POOL_PREFIX}.PureCLIP.crosslink.bed \
#  -b 10 -g ${GENOME_SIZE} | bedtools shuffle -i stdin \
#    -g ${GENOME_SIZE} -noOverlapping -seed 1000 -chrom | \
#      bedtools getfasta -fi ${FASTA} -bed stdin -s \
#        -name+ > ${POOL_PREFIX}.PureCLIP.negative.fa

module load MEME/5.1.0-OpenMPI.3.0.0
#dreme -p ${POOL_PREFIX}.PureCLIP.postive.fa \
#  -norc -k 6 -m 10 -s 1000 -png -oc ${POOL_PREFIX}_DREME -e 0.1 \
#  > ${POOL_PREFIX}.dreme.log 2>&1
#
#export WINEXTRACT=`which winextract`
#compute_CLmotif_scores.sh ${FASTA} ${MERGE_PREFIX}.bam \
#  ${POOL_PREFIX}_DREME/dreme.xml ${POOL_PREFIX}_DREME/dreme.txt \
#  ${MERGE_PREFIX}.fimo_clmotif_occurences.bed

pureclip -i ${MERGE_PREFIX}.bam \
  -bai ${MERGE_PREFIX}.bam.bai \
  -g ${FASTA} -iv 'chr1;chr2;chr3;' -nt ${THREAD} \
  -o ${POOL_PREFIX}.PureCLIP.crosslink.cov_CLmotifs.bed \
  -or ${POOL_PREFIX}.PureCLIP.bind.cov_CLmotifs.bed \
  -p ${POOL_PREFIX}.PureCLIP.parameter.cov_CLmotifs.txt \
  -nim 4 -fis ${MERGE_PREFIX}.fimo_clmotif_occurences.bed \
  > ${MERGE_PREFIX}.pureCLIP.cov_CLmotifs.log 2>&1

#
#bedtools slop -i ${POOL_PREFIX}.PureCLIP.crosslink.bed \
#  -b 5 -g ${GENOME_SIZE} > ${POOL_PREFIX}.temp.bed
#
#scanMotif.py -input ${POOL_PREFIX}.temp.bed -format bed6 \
#  -fasta ${FASTA} -motif RRACH -tag 3 \
#  -output ${POOL_PREFIX}.PureCLIP.crosslink.RRACH.all.bed
#
#bedtools shift -i ${POOL_PREFIX}.PureCLIP.crosslink.RRACH.all.bed \
#  -g ${GENOME_SIZE} -p 1 -m -1 > ${POOL_PREFIX}.temp.bed
#
#ctk_C2T_mutation_filter.sh "${POOL_PREFIX}.temp.bed" \
#  "${POOL_PREFIX}.pooled.mutation.txt" \
#  "${POOL_PREFIX}.PureCLIP.crosslink.RRACH.CT.bed"
#
#bedtools shift -i ${POOL_PREFIX}.PureCLIP.crosslink.RRACH.CT.bed \
#  -g ${GENOME_SIZE} -m 1 -p -1 > ${POOL_PREFIX}.PureCLIP.crosslink.RRACH.bed
#
#rm -f ${POOL_PREFIX}.temp.bed \
#  ${POOL_PREFIX}.PureCLIP.crosslink.RRACH.CT.bed \
#  ${POOL_PREFIX}.PureCLIP.crosslink.RRACH.all.bed
