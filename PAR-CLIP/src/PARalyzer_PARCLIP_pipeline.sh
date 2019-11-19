#!/bin/bash
#SBATCH --job-name=PARalyzer_PARCLIP_pipeline    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 10                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=100G                      # Amount of memory in GB
#SBATCH --time=72:10:00               # Time limit hrs:min:sec
#SBATCH --output=PARalyzer_PARCLIP_pipeline.log               # Time limit hrs:min:sec

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
# ref:https://ohlerlab.mdc-berlin.de/files/duke/PARalyzer/README.txt

function showHelp {
  echo -ne "usage: sbatch PARalyzer_PARCLIP_pipeline.sh -n <thread_num> -o <log> --mem <200G> "
  echo -ne "PARalyzer_PARCLIP_pipeline.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -b | --barcode-length: barcode length <int>
    -e | --exp-prefix: experiment prefix string <str>
    -f | --fasta: genome fasta file <str>
    -g | --gsize: genome size file <str>
    -i | --input: input fastq directory (cutadapt) <str>
    -o | --output: output result directory <str>
    -p | --pool-prefix: pooled experiment prefix <str>
    -t | --thread: # of cpus <int>
    --2bit: Genome .2bit file (UCSC) <str>
    --gtf: genome annotation gtf <str>
    --full-bed: all transcripts annotation in bed12 <str>
    --longest-bed: mRNA longest annotation in bed12 <str>
    --max-mismatch: --outFilterMismatchNoverLmax in STAR (0.1) <float>
    --mem: # of memory used for PARalyzer (20G) <int>
    --CLcount: # of ConversionLocationCount (1) <int>
    --CEcount: # of ConversionEventCount (1) <int>
    --Rcount: # of ReadCount (1) <int>
    --nfrom: modified ribonucleotide (A|T|C|T) (T) <str>
    --nto: converted ribonucleotide (A|T|C|T) (C) <str>
    --repeat-bed: repeat bed used for filtering (eg.t/rRNA) <str>
    --downstream: search MOTIF at #bp downstream of transition site <str>
    --upstream: search MOTIF at #bp upstream of transition site <str>
    --skip-mapping: skip reads mapping step <bool>
    --skip-calling: skip crosslink sites calling step <bool>
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

TEMP=`getopt -o hb:e:g:i:o:p:t: --long help,skip-mapping,skip-calling \
  --long bowtie,STAR,PCR,keep-tmp-fastq \
  --long input:,thread:,output:,exp-prefix:,pool-prefix:,index: \
  --long longest-bed:,full-bed:,gtf:,repeat-bed:,barcode-length:,fasta: \
  --long nfrom:,nto:,mem:,motif:,mtag:,downstream:,upstream: \
  --long 2bit:,CLcount:,CEcount:,Rcount: \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREAD=10
MEMORY="20G"
BARCODE_LEN=9
OUTPUT_DIR=
INPUT_DIR=
EXP_PREFIX=
POOL_PREFIX=
GENOME_INDEX=
FASTA=
GENOME_2BIT=
GTF=
GENOME_SIZE=
LONGEST_BED=
FULL_BED=
REPEAT_BED=
MAX_MISMATCH=0.1
CL_COUNT=1
CE_COUNT=1
R_COUNT=1
N_FROM="T"
N_TO="C"
MOTIF="RRACH"
MOTIF_TAG=3
DOWNSTREAM=10
UPSTREAM=10
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
    --2bit ) GENOME_2BIT="$2"; shift 2 ;;
    --CLcount ) CL_COUNT="$2"; shift 2 ;;
    --CEcount ) CE_COUNT="$2"; shift 2 ;;
    --Rcount ) R_COUNT="$2"; shift 2 ;;
    --gtf ) GTF="$2"; shift 2 ;;
    --full-bed ) FULL_BED="$2"; shift 2 ;;
    --index ) GENOME_INDEX="$2"; shift 2 ;;
    --longest-bed ) LONGEST_BED="$2"; shift 2 ;;
    --max-mismatch ) MAX_MISMATCH="$2"; shift 2 ;;
    --mem ) MEMORY="$2"; shift 2 ;;
    --motif ) MOTIF="$2"; shift 2 ;;
    --mtag ) MOTIF_TAG="$2"; shift 2 ;;
    --nfrom ) N_FROM="$2"; shift 2 ;;
    --nto ) N_TO="$2"; shift 2 ;;
    --repeat-bed ) REPEAT_BED="$2"; shift 2 ;;
    --keep-tmp-fastq ) KEEP_TMP_FASTQ=true; shift ;;
    --skip-mapping ) SKIP_MAPPING=true; shift ;;
    --skip-calling ) SKIP_CALLING=true; shift ;;
    --downstream ) DOWNSTREAM="$2"; shift 2 ;;
    --upstream ) UPSTREAM="$2"; shift 2 ;;
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

if [ -z $GENOME_2BIT ]; then
  echo "--2bit: Wrong! genome.2bit NOT found!"
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

if [[ $N_FROM == $N_TO ]]; then
  echo "--nfrom and --nto: Wrong! They are should not be the same!"
  exit 2
fi

echo "Running pipleline with following parameters:"
echo "THREAD=$THREAD"
echo "MEMORY=$MEMORY"
echo "BARCODE_LEN=$BARCODE_LEN"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo "INPUT_DIR=$INPUT_DIR"
echo "EXP_PREFIX=$EXP_PREFIX"
echo "POOL_PREFIX=$POOL_PREFIX"
echo "GENOME_INDEX=$GENOME_INDEX"
echo "FASTA=$FASTA"
echo "GENOME_2BIT=$GENOME_2BIT"
echo "GTF=$GTF"
echo "GENOME_SIZE=$GENOME_SIZE"
echo "LONGEST_BED=$LONGEST_BED"
echo "FULL_BED=$FULL_BED"
echo "REPEAT_BED=$REPEAT_BED"
echo "MAX_MISMATCH=$MAX_MISMATCH"
echo "CL_COUNT=$CL_COUNT"
echo "CE_COUNT=$CE_COUNT"
echo "R_COUNT=$R_COUNT"
echo "N_FROM=$N_FROM"
echo "N_TO=$N_TO"
echo "MOTIF=$MOTIF"
echo "MOTIF_TAG=$MOTIF_TAG"
echo "DOWNSTREAM=$DOWNSTREAM"
echo "UPSTREAM=$UPSTREAM"
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
RESULT_DIR="$OUTPUT_DIR/cluster"
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

## Read preprocessing and mapping
if [[ ! -d $MAP_DIR  ]]; then
  mkdir -p $MAP_DIR
fi

if [[ ! -d $FASTQC_DIR ]]; then
  mkdir $FASTQC_DIR
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
    {
      cd $MAP_DIR
      PREFIX=${i%%.fastq}
      PREFIX=${PREFIX##*/} # *.trim
      MAP_PREFIX=${i%%.trim.fastq}
      MAP_PREFIX=${MAP_PREFIX##*/}

      if $STAR_FLAG; then
        if $PCR_FLAG; then
          echo "${MAP_PREFIX}: Activate --PCR! PCR duplicates removing..."
          ## collapse PCR duplicates
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
        bowtie -p ${MAP_THREAD} -t -v 2 -m 10 --best --strata --sam $GENOME_INDEX \
          ${input} ${MAP_PREFIX}.aligned.sam > ${MAP_PREFIX}.aligned.log 2>&1
        ## sam -> bam, sort bam
        samtools view -@ ${MAP_THREAD} -h -bS -F 4 \
          ${MAP_PREFIX}.aligned.sam -o ${MAP_PREFIX}.aligned.unsort.bam
        samtools sort -@ ${MAP_THREAD} -m 2G -O bam \
          -o ${MAP_PREFIX}.aligned.bam ${MAP_PREFIX}.aligned.unsort.bam
        ## index bam
        samtools index -@ ${MAP_THREAD} ${MAP_PREFIX}.aligned.bam
        ## delete sam
        rm -f ${MAP_PREFIX}.aligned.sam ${MAP_PREFIX}.aligned.unsort.bam
        ## fastqc check mapped bam
        fastqc -o $FASTQC_DIR ${POOL_PREFIX}.aligned.bam \
          > ${POOL_PREFIX}.aligned.fastqc.log 2>&1
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
fi

## convert bam to sam
echo "Converting bam to sam..."
cd $MAP_DIR

if [[ ! -d $RESULT_DIR  ]]; then
  mkdir -p $RESULT_DIR
fi

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
  {
    PREFIX=${i%%.fastq}
    PREFIX=${PREFIX##*/} # *.trim
    MAP_PREFIX=${i%%.trim.fastq}
    MAP_PREFIX=${MAP_PREFIX##*/}
    samtools view -@ ${MAP_THREAD} -h -O SAM \
      -o ${MAP_PREFIX}.aligned.sam ${MAP_PREFIX}.aligned.bam
    ### link sam to $RESULT_DIR
    ln -sf $MAP_DIR/${MAP_PREFIX}.aligned.sam $RESULT_DIR
    echo >&9
  } &
done
wait
exec 9>&-
exec 9<&-

## PARalyzer start

### functions
function setupini () {
  INI="BANDWIDTH=3\n"
  INI+="CONVERSION=${N_FROM}>${N_TO}\n"
  INI+="MINIMUM_READ_COUNT_PER_GROUP=10\n"
  INI+="MINIMUM_READ_COUNT_PER_CLUSTER=5\n"
  INI+="MINIMUM_READ_COUNT_FOR_KDE=5\n"
  INI+="MINIMUM_CLUSTER_SIZE=8\n"
  INI+="MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER=1\n"
  INI+="MINIMUM_CONVERSION_COUNT_FOR_CLUSTER=1\n"
  INI+="MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION=1\n"
  INI+="MINIMUM_READ_LENGTH=13\n"
  INI+="MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES=5\n"
  INI+="\n"
  INI+="EXTEND_BY_READ\n"
  INI+="\n"
  INI+="${INPUT_ALIGNMENTS}"
  INI+="\n"
  INI+="GENOME_2BIT_FILE=${GENOME_2BIT}\n"
  INI+="\n"
  INI+="${FILTER}"
  INI+="OUTPUT_CLUSTERS_FILE=${POOL_PREFIX}.cluster.csv\n"
  INI+="OUTPUT_GROUPS_FILE=${POOL_PREFIX}.groups.csv\n"
  INI+="OUTPUT_DISTRIBUTIONS_FILE=${POOL_PREFIX}_distribution.csv\n"
  printf "%b" "${INI}" > ${INI_FILE}
}

## calling crosslink sites
if $SKIP_CALLING; then
  echo "Skip calling crosslink sites."
else
  if [[ ! -d $RESULT_DIR  ]]; then
    mkdir -p $RESULT_DIR
  fi
  cd $RESULT_DIR

  ### INPUT_ALIGNMENTS
  INPUT_ALIGNMENTS=""
  for SAM_FILE in `find ./ -maxdepth 1 -type l -name "${EXP_PREFIX}*.aligned.sam" | sort`;
  do
    INPUT_ALIGNMENTS+="SAM_FILE=${SAM_FILE}\n"
  done
  
  ### setting INI file
  INI_FILE="${POOL_PREFIX}.PARalyzer.ini"
  if [[ ! -z $REPEAT_BED ]]; then
    flag=${REPEAT_BED##*/}
    flag=${flag%%.*}
    FILTER="FILTER_FILE=${REPEAT_BED}=${flag}\n"
  else
    FILTER=""
  fi

  if [[ ! -f $GENOME_2BIT ]]; then
    echo "No genome.2bit file found!"
    echo "Generating genome.2bit file..."
    faToTwoBit ${FASTA} ${GENOME_2BIT}
    echo "Genome.2bit is ready."
  fi
  
  ## generating INI file
  setupini
  ## start PARalyzer
  echo "PARalyzer: crosslink sites calling..."
  PARalyzer ${MEMORY} ${INI_FILE} > ${POOL_PREFIX}.PARalyzer.log 2>&1
  echo "PARalyzer calling crosslink sites done."
fi

## Pooling CT and Truncation m6A sites
if [[ ! -d $FINAL_DIR  ]]; then
  mkdir -p $FINAL_DIR
fi
cd $FINAL_DIR

### link cluster.csv and groups.csv
ln -sf ${RESULT_DIR}/${POOL_PREFIX}.cluster.csv ./
ln -sf ${RESULT_DIR}/${POOL_PREFIX}.groups.csv ./

echo "Converting coordinates to bed-format..."

if [[ -z $REPEAT_BED ]]; then
  repeat=1
else
  repeat=0
fi

awk -v clc="$CL_COUNT" -v cec="$CE_COUNT" -v rc="$R_COUNT" \
  -v rp="$repeat" 'BEGIN { FS=","; OFS="\t";}
{
  if (rp==1) {
    if(FNR>1 && $13=="NA") {
      if ( $10>=clc && $11>=cec && $7>=rc ) {
        $3=$3-1; if($3<0){$3=0};
        print $1,$3,$4,$5,$7,$2,$6,$8,$9,$10,$11,$12;
      }
    }
  }else{
    if (FNR>1) {
      if ( $10>=clc && $11>=cec && $7>=rc ) {
        $3=$3-1; if($3<0){$3=0};
        print $1,$3,$4,$5,$7,$2,$6,$8,$9,$10,$11,$12;
      }
    }
  }
}' ${POOL_PREFIX}.cluster.csv | \
  sort -t $'\t' -k 1,1 -k 2,2n > ${POOL_PREFIX}.cluster.bed

awk -v clc="$CL_COUNT" -v cec="$CE_COUNT" -v rc="$R_COUNT" \
  -v rp="$repeat" 'BEGIN { FS=","; OFS="\t";}
{
  if (rp==1) {
    if(FNR>1 && $13=="NA") {
      if ( $8>=clc && $9>=cec && $7>=rc ) {
        $3=$3-1; if($3<0){$3=0};
        print $1,$3,$4,$5,$7,$2,$6,$8,$9;
      }
    }
  }else{
    if (FNR>1) {
      if ( $8>=clc && $9>=cec && $7>=rc ) {
        $3=$3-1; if($3<0){$3=0};
        print $1,$3,$4,$5,$7,$2,$6,$8,$9;
      }
    }
  }
}' ${POOL_PREFIX}.groups.csv | \
  sort -t $'\t' -k 1,1 -k 2,2n > ${POOL_PREFIX}.groups.bed

echo "csv conversion done."

if [[ -z $MOTIF ]]; then
  ## scan ${MOTIF} motif in peaks
  echo "Scanning ${MOTIF} motif in peaks..."
  peaks=( ${POOL_PREFIX}.cluster.bed ${POOL_PREFIX}.groups.bed)
  for i in "${peaks[@]}";
  do
    prefix=${i%%.bed}
    ## keep peaks with length >= ${SPAN_WIDTH}
    awk -v span="$SPAN_WIDTH" 'BEGIN{FS="\t";OFS="\t";}{
      if(($3-$2)>=span){print}}' $i > $i.keep.tmp
    awk -v span="$SPAN_WIDTH" 'BEGIN{FS="\t";OFS="\t";}{
      if(($3-$2)<span){
        $2 = int(($3+$2)/2);
        $3 = start + 1;
        print $0;
      }
    }' $i | bedtools slop -i ${i} -l ${UPSTREAM} \
      -r ${DOWNSTREAM} -s -g ${GENOME_SIZE} > $i.extend.tmp
    cat $i.keep.tmp $i.extend.tmp | sort -k1,1 -k2,2n \
      > $i.tmp
  
    scanMotif.py -input $i.tmp -format bed6 \
      -fasta ${FASTA} -motif ${MOTIF} -tag ${MOTIF_TAG} \
      -output ${prefix}.${MOTIF}.bed
  done
  rm -f *.tmp
fi

## annotate beds
echo "Annotating beds..."
if [ ! -z $LONGEST_BED ]; then
  for i in `find ./ -type f -name "${POOL_PREFIX}*.bed"`;
  do
    temp="${POOL_PREFIX}.tmp"
    cut -f 1-6 $i > $temp
    PREFIX=${i%%.bed}
    bedBinDistribution.pl -input $temp -bed12 $LONGEST_BED \
      --type count -o ${PREFIX}.count.bin
    bedBinDistribution.pl -input $temp -bed12 $LONGEST_BED \
      -o ${PREFIX}.percentage.bin
    paste ${PREFIX}.count.bin ${PREFIX}.percentage.bin | cut -f 1,2,3,6 > ${PREFIX}.bin
    sed -i '1i region\tbin\tCount\tPercentage' ${PREFIX}.bin
    rm -f ${PREFIX}.count.bin ${PREFIX}.percentage.bin
  done
fi

if [ ! -z $FULL_BED ]; then
  ## get mRNA annotation bed12
  cat $FULL_BED | awk '/protein_coding.+protein_coding\t/' > mRNA.annotation.bed12.tmp
  for i in `find ./ -type f -name "${POOL_PREFIX}*.bed"`;
  do
    temp="${POOL_PREFIX}.tmp"
    cut -f 1-6 $i > $temp
    PREFIX=${i%%.bed}
    ### gene type
    geneDistribution.pl -strand --input $temp \
      -bed12 $FULL_BED -o ${PREFIX}.gene
    sed -i '1i geneType\tpeakNumber' ${PREFIX}.gene
    ### gene region
    regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' \
      --input $temp \
      -bed12 mRNA.annotation.bed12.tmp -o ${PREFIX}.region
    sed -i '1i region\tpeakNumber\tenrichment' ${PREFIX}.region
  done
  rm -f mRNA.annotation.bed12.tmp
fi
rm -f *.tmp

echo "beds annotation done."
