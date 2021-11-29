#!/bin/bash
#SBATCH --job-name=miCLIP_rbsSeeker_pipeline    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 10                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=100G                      # Amount of memory in GB
#SBATCH --time=120:10:00               # Time limit hrs:min:sec
#SBATCH --output=miCLIP_rbsSeeker_pipeline.log               # Time limit hrs:min:sec

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.

function showHelp {
  echo -ne "usage: sbatch miCLIP_rbsSeeker_pipeline.sh -n <thread_num> -o <log> --mem <200G> "
  echo -ne "miCLIP_rbsSeeker_pipeline.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -b | --barcode-length: barcode length <int>
    -e | --exp-prefix: experiment prefix string <str>
    -f | --fasta: genome fasta file <str>
    -g | --gsize: genome size file <str>
    -i | --input: input fastq directory (exp_prefix*.trim.fastq) <str>
    -o | --output: output result directory <str>
    -p | --pool-prefix: pooled experiment prefix <str>
    -t | --thread: # of cpus <int>
    --gtf: genome annotation gtf <str>
    --full-bed: all transcripts annotation in bed12 <str>
    --max-mismatch: --outFilterMismatchNoverLmax in STAR <float>
    --repeat-bed: repeat bed used for filtering (eg.t/rRNA) <str>
    --skip-mapping: skip reads mapping step <bool>
    --skip-calling: skip m6A sites calling step <bool>
    --bowtie: map reads with bowtie aligner <bool>
    --PCR: Collapse reads before mapping (CTK:fastq2collapse.pl) <bool>
    --keep-tmp-fastq: keep temporary fastqs <bool>
    --STAR: map reads with STAR aligner <bool>
    --index: genome index <str>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o hb:e:g:i:o:p:t:, --long help,skip-mapping,skip-calling, \
  --long bowtie,STAR,PCR,keep-tmp-fastq, \
  --long input:,thread:,output:,exp-prefix:,pool-prefix:,index:, \
  --long full-bed:,gtf:,repeat-bed:,barcode-length:,fasta:, \
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
FULL_BED=
REPEAT_BED=
MAX_MISMATCH=0.1
BOWTIE_FLAG=false
PCR_FLAG=false
KEEP_TMP_FASTQ=false
STAR_FLAG=false
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
    --max-mismatch ) MAX_MISMATCH="$2"; shift 2 ;;
    --repeat-bed ) REPEAT_BED="$2"; shift 2 ;;
    --keep-tmp-fastq ) KEEP_TMP_FASTQ=true; shift ;;
    --skip-mapping ) SKIP_MAPPING=true; shift ;;
    --skip-calling ) SKIP_CALLING=true; shift ;;
    --bowtie ) BOWTIE_FLAG=true; shift ;;
    --PCR ) PCR_FLAG=true; shift ;;
    --STAR ) STAR_FLAG=true; shift ;;
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
  echo "-o|--output: Wrong! Output directory NOT found!"
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

if [ -z $STAR_INDEX ]; then
  if [ -z $FASTA ]; then
    echo "-f|--fasta: Wrong! Genome fasta NOT found!"
    exit 2
  fi
  if [ -z $GTF ]; then
    echo "--gtf: Wrong! Annotation gtf NOT found!"
    exit 2
  fi
fi

echo "Running pipleline with following parameters:"
echo "THREAD=$THREAD"
echo "BARCODE_LEN=$BARCODE_LEN"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo "INPUT_DIR=$INPUT_DIR"
echo "EXP_PREFIX=$EXP_PREFIX"
echo "POOL_PREFIX=$POOL_PREFIX"
echo "GENOME_INDEX=$GENOME_INDEX"
echo "FASTA=$FASTA"
echo "GTF=$GTF"
echo "GENOME_SIZE=$GENOME_SIZE"
echo "FULL_BED=$FULL_BED"
echo "REPEAT_BED=$REPEAT_BED"
echo "MAX_MISMATCH=$MAX_MISMATCH"
echo "BOWTIE_FLAG=$BOWTIE_FLAG"
echo "PCR_FLAG=$PCR_FLAG"
echo "KEEP_TMP_FASTQ=$KEEP_TMP_FASTQ"
echo "STAR_FLAG=$STAR_FLAG"
echo "SKIP_MAPPING=$SKIP_MAPPING"
echo "SKIP_CALLING=$SKIP_CALLING"
echo ""

# basic variables
MAP_DIR="$OUTPUT_DIR/mapping"
FASTQC_DIR="$OUTPUT_DIR/mapping/fastQC"
RESULT_DIR="$OUTPUT_DIR/result"
FINAL_DIR="$OUTPUT_DIR/final"

REP_NUM=`find $INPUT_DIR -type f -name "${EXP_PREFIX}*.trim.fastq" | wc -l`
MAP_THREAD=$((THREAD / REP_NUM))

if (( $MAP_THREAD == 0 )); then
  MAP_THREAD=1
fi

if [[ ! -d $OUTPUT_DIR ]]; then
  mkdir -p $OUTPUT_DIR
fi
## try to build index
if ! $SKIP_MAPPING; then
  #statements
  if $STAR_FLAG; then
    echo "Mapping with STAR."
    if [ -z $GENOME_INDEX ]; then
      GENOME_INDEX="$OUTPUT_DIR/STAR_index/${POOL_PREFIX}"
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
  if $BOWTIE_FLAG; then
    version=`bowtie --version`
    echo "Mapping with bowtie: ${version}"
    if [[ ! -f ${GENOME_INDEX}.1.ebwt ]]; then
      echo "NO bowtie index detected!"
      exit 2
    fi
  fi
fi
## genome

## Read preprocessing and mapping

if [[ ! -d $MAP_DIR  ]]; then
    mkdir -p $MAP_DIR
fi

if $SKIP_MAPPING ; then
  echo "Skip mapping."
else
  echo "Mapping start..."
  ### parallel start
  tmpfifo="/tmp/$RANDOM.fifo"
  mkfifo tmpfifo
  exec 9<>tmpfifo
  rm -rf tmpfifo
  
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

      if $STAR_FLAG; then
        if $PCR_FLAG; then
          echo "${MAP_PREFIX}: Activate --PCR! PCR duplicates removing..."
          ## remove PCR duplicates
          fastq2collapse.pl $i - | gzip -c > ${PREFIX}.c.fastq.gz
          input=${PREFIX}.c.fastq.gz
        else
          input=$i
        fi
        if (( $BARCODE_LEN > 0 )); then
          ## strip barcode
          stripBarcode.pl -format fastq -len ${BARCODE_LEN} ${input} - | \
            gzip -c > ${PREFIX}.bc.fastq.gz
          input=${PREFIX}.bc.fastq.gz
        fi
        echo "${MAP_PREFIX}: Mapping reads with STAR..."
        ## mapping, reduce mismatch rates
        STAR --outSAMtype BAM SortedByCoordinate --runThreadN ${MAP_THREAD} \
          --genomeDir ${GENOME_INDEX} --readFilesIn ${input} \
          --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --alignEndsType EndToEnd \
          --outFilterType Normal --outFilterMultimapScoreRange 0 \
          --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
          --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
          --outFilterMismatchNoverLmax ${MAX_MISMATCH}  --outFilterScoreMin 0 \
          --outFilterScoreMinOverLread 0 --outFilterMatchNmin 15  \
          --outFilterMatchNminOverLread 0 --alignIntronMin 1 --alignIntronMax 1 \
          --alignMatesGapMax 1500  --seedSearchStartLmax 15 \
          --seedSearchStartLmaxOverLread 1 --seedSearchLmax 0 --seedMultimapNmax 20000 \
          --seedPerReadNmax 1000 --seedPerWindowNmax 100 --seedNoneLociPerWindow 20 \
          --outSAMmode Full --outSAMattributes All --outSAMunmapped None \
          --outSAMorder Paired --outSAMprimaryFlag AllBestScore \
          --outSAMreadID Standard --outReadsUnmapped None --readFilesCommand zcat \
          --outFileNamePrefix ${MAP_PREFIX}. --outTmpDir ${MAP_PREFIX}_STARtmp
        rm -rf ${MAP_PREFIX}_STARtmp
        # index
        #mv ${MAP_PREFIX}.Aligned.sortedByCoord.out.bam ${MAP_PREFIX}.aligned.bam
        #samtools index -@ ${MAP_THREAD} ${MAP_PREFIX}.aligned.bam
        samtools view -@ ${MAP_THREAD} -hb ${MAP_PREFIX}.Aligned.sortedByCoord.out.bam \
          -o ${MAP_PREFIX}.aligned.bam
        samtools index -@ ${MAP_THREAD} ${MAP_PREFIX}.aligned.bam
        rm -f ${MAP_PREFIX}.Aligned.sortedByCoord.out.bam ${MAP_PREFIX}.Aligned.sortedByCoord.out.bam.sai
      else
        if $PCR_FLAG; then
          echo "${MAP_PREFIX}: Activate --PCR! PCR duplicates removing..."
          ## remove PCR duplicates
          fastq2collapse.pl $i ${PREFIX}.c.fastq
          input=${PREFIX}.c.fastq
        else
          input=$i
        fi
        if (( $BARCODE_LEN > 0 )); then
          ## strip barcode
          stripBarcode.pl -format fastq -len ${BARCODE_LEN} ${input} ${PREFIX}.bc.fastq
          input=${PREFIX}.bc.fastq
        fi
        ## bowtie alignment
        echo "${MAP_PREFIX}: Mapping reads with bowtie..."
        bowtie -p ${MAP_THREAD} -t -v 2 -m 20 --best --strata --sam $GENOME_INDEX \
          ${input} ${MAP_PREFIX}.aligned.sam > ${MAP_PREFIX}.aligned.log 2>&1
        samtools view -@ ${MAP_THREAD} -h -bS -F 4 \
          ${MAP_PREFIX}.aligned.sam -o ${MAP_PREFIX}.aligned.unsort.bam
        samtools sort -@ ${MAP_THREAD} -m 2G -O bam \
          -o ${MAP_PREFIX}.aligned.bam ${MAP_PREFIX}.aligned.unsort.bam
        rm -f ${MAP_PREFIX}.aligned.sam ${MAP_PREFIX}.aligned.unsort.bam
        samtools index -@ ${MAP_THREAD} ${MAP_PREFIX}.aligned.bam
        echo "${MAP_PREFIX}: Reads mapping done."
      fi
      echo >&9
    } &
  done
  wait
  exec 9>&-
  exec 9<&-
  # parallel end
  echo "Reads mapping of all replicates are done."
  cd $MAP_DIR
  #delete tmp fastqs
  if $KEEP_TMP_FASTQ; then
    echo "Keep tmp fastqs."
  else
    echo "Deleting tmp fastqs..."
    find ./ -maxdepth 1 -type f -name "${EXP_PREFIX}*.fastq" | xargs -I {} rm -f {}
    find ./ -maxdepth 1 -type f -name "${EXP_PREFIX}*.fastq.gz" | xargs -I {} rm -f {}
  fi
  ### pooling reads
  echo "Polling reads..."
  DUPLRM_BAMS=""
  for i in `find ./ -maxdepth 1 -type f -name "${EXP_PREFIX}*.aligned.bam" | sort`;
  do
    DUPLRM_BAMS=`echo "${DUPLRM_BAMS}"" ""${i}"`
  done
  
  samtools merge -@ ${MAP_THREAD} -f ${POOL_PREFIX}.aligned.pooled.bam ${DUPLRM_BAMS}
  samtools index -@ ${MAP_THREAD} ${POOL_PREFIX}.aligned.pooled.bam
  
  ### Quality control
  if [[ ! -d $FASTQC_DIR ]]; then
      mkdir $FASTQC_DIR
  fi
  fastqc -o $FASTQC_DIR ${POOL_PREFIX}.aligned.pooled.bam \
    > ${POOL_PREFIX}.aligned.pooled.fastqc.log 2>&1
  echo "Reads pooled."
fi

## calling crosslink sites
if $SKIP_CALLING; then
  echo "Skip calling m6A sites."
else
  if [[ ! -d $RESULT_DIR  ]]; then
      mkdir -p $RESULT_DIR
  fi
  cd $RESULT_DIR
  
  ln -sf ${MAP_DIR}/${POOL_PREFIX}.aligned.pooled.bam.bai ./
  ln -sf ${MAP_DIR}/${POOL_PREFIX}.aligned.pooled.bam ./
  
  INPUT_BAM="${POOL_PREFIX}.aligned.pooled.bam"
  
  # rbsSeeker m6A sites calling
  echo "rbsSeeker m6A sites calling..."
  rbsSeeker -T CT -L 20 -t 129600000 -n 1 -H 3 -d 1 \
    -p 1 -q 1 -o "$RESULT_DIR" -P "$POOL_PREFIX" \
    --fa "${FASTA}" --fai "${FASTA}.fai" \
    --bam "$INPUT_BAM" \
    > ${POOL_PREFIX}.pooled.rbsSeeker.log 2>&1
  
  echo "rbsSeeker calling m6A sites done."
fi

## Pooling CT and Truncation m6A sites
if [[ ! -d $FINAL_DIR  ]]; then
  mkdir -p $FINAL_DIR
fi

cd $FINAL_DIR

### link CT and Truncation beds
ln -sf ${RESULT_DIR}/${POOL_PREFIX}_rbsSeeker_CT.bed ./${POOL_PREFIX}_rbsSeeker_CIMS.bed
ln -sf ${RESULT_DIR}/${POOL_PREFIX}_rbsSeeker_Truncation.bed ./${POOL_PREFIX}_rbsSeeker_CITS.bed

echo "Pooling CT and Truncation m6A sites..."

for i in `find ./ -type l -name "${POOL_PREFIX}*.bed" | grep -E "rbsSeeker_(CIMS|CITS)"`;
do
  outputName="${i//_rbsSeeker_/.rbsSeeker.}"
  awk 'BEGIN{OFS="\t";FS="\t";}
  {
    if(FNR>1){
      if($8 == 7){
        seq = substr($7, 8, 5);
        print $1,$2,$3,$4,int($12),$6,$10,$14,seq;
      }
    }
  }' $i | sort -k1,1 -k2,2n | \
  bedtools shift -i stdin -g ${GENOME_SIZE} -m 1 -p -1 > ${POOL_PREFIX}.tmp.bed
  ## filter *.bed file and get coordinates corresponding to "A"
  getBedBySeq.py -fasta $FASTA -bed ${POOL_PREFIX}.tmp.bed -seq "A" -output $outputName
  unlink $i
  rm -f ${POOL_PREFIX}.tmp.bed
done

prefixArr=( "${POOL_PREFIX}" )
for i in "${prefixArr[@]}"
do
  cat ${i}.rbsSeeker.CITS.bed ${i}.rbsSeeker.CIMS.bed | \
    awk 'BEGIN{OFS="\t";FS="\t";}
    {
      key = $1"\t"$2"\t"$3":"$6;
      if (key in arrayA) {
        arrayA[key] = arrayA[key]"|"$4
        arrayB[key] = arrayB[key] + int($5)
      }else{
        arrayA[key] = $4
        arrayB[key] = int($5)
        arrayC[key] = $9
      }
    }
    END{
      for (key in arrayA) {
        split(key,splitArr,":");
        pos = splitArr[1];
        strand = splitArr[2];
        print pos, arrayC[key]"="arrayA[key], arrayB[key], strand;
      }
    }' | sort -k1,1 -k2,2n > ${i}.rbsSeeker.combine.bed
done

echo "Pooling CT and Truncation m6A sites done."

## annotate beds
echo "Annotating beds..."
BED_FILES=""
for bed in `find ./ -type f -name "${POOL_PREFIX}*.bed"`;
do
  BED_FILES="$BED_FILES $bed"
done

metagene.py --anno $FULL_BED -i $BED_FILES \
  --method center -b constant --smooth move --output ${POOL_PREFIX}.metagene.txt

#if [ ! -z $FULL_BED ]; then
#  ## get mRNA annotation bed12
#  cat $FULL_BED | awk '/protein_coding.+protein_coding\t/' > mRNA.annotation.bed12.tmp
#  for i in `find ./ -type f -name "${POOL_PREFIX}*.bed"`;
#  do
#    temp="${POOL_PREFIX}.tmp"
#    cut -f 1-6 $i > $temp
#    PREFIX=${i%%.bed}
#    ### gene type
#    geneDistribution.pl -strand --input $temp \
#      -bed12 $FULL_BED -o ${PREFIX}.gene
#    sed -i '1i geneType\tpeakNumber' ${PREFIX}.gene
#    ### gene region
#    regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' \
#      --input $temp \
#      -bed12 mRNA.annotation.bed12.tmp -o ${PREFIX}.region
#    sed -i '1i region\tpeakNumber\tenrichment' ${PREFIX}.region
#  done
#  rm -f mRNA.annotation.bed12.tmp
#fi
#rm -f *.tmp

echo "beds annotation done."
