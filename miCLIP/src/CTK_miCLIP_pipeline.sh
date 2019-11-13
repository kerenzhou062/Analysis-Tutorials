#!/bin/bash
#SBATCH --job-name=CTK_miCLIP_pipeline    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail  
#SBATCH -n 10                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=100G                      # Amount of memory in GB
#SBATCH --time=72:10:00               # Time limit hrs:min:sec
#SBATCH --output=CTK_miCLIP_pipeline.log   # Standard output and error log

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.

function showHelp {
  echo -ne "usage: sbatch CTK_miCLIP_pipeline.sh -n <thread_num> -o <log> --mem <200G> "
  echo -ne "CTK_miCLIP_pipeline.sh <options>\n"
  echo -e "options:
    -h | --help: show help infomation <bool>
    -b | --barcode-length: barcode length <int>
    -e | --exp-prefix: experiment prefix string <str>
    -f | --fasta: genome fasta file <str>
    -g | --gsize: genome size file <str>
    -i | --input: input fastq directory (cutadapt) <str>
    -m | --min-length: minimun length of mapped reads to parse (parseAlignment) <str>
    -o | --output: output result directory <str>
    -p | --pool-prefix: pooled experiment prefix <str>
    -q | --quality: quality cutoff of mapped reads to parse (parseAlignment) <str>
    -t | --thread: # of cpus <int>
    --dbkey: CTK dbkey (hg38|hg19|mm10) <str>
    --full-bed: all transcripts annotation in bed12 <str>
    --index: genome index <str>
    --longest-bed: mRNA longest annotation in bed12 <str>
    --repeat-bed: repeat bed used for filtering (eg.t/rRNA) <str>
    --skip-mapping: skip reads mapping step <bool>
    --skip-pooling: skip pooling step <bool>
    --skip-calling: skip m6A sites calling step <bool>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]
then
  showHelp
  exit 2
fi

TEMP=`getopt -o hb:e:g:i:m:o:p:q:t: --long help,skip-mapping,skip-pooling,skip-calling \
  --long input:,thread:,output:,exp-prefix:,pool-prefix:,index:,min-length:,dbkey: \
  --long longest-bed:,full-bed:,repeat-bed:,barcode-length:,fasta:,quality: \
  -- "$@"`


if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREAD=10
BARCODE_LEN=9
MIN_LENGTH=18
OUTPUT_DIR=
INPUT_DIR=
EXP_PREFIX=
POOL_PREFIX=
BWA_INDEX=
DB_KEY=
QUALITY=20
FASTA=
GENOME_SIZE=
LONGEST_BED=
FULL_BED=
REPEAT_BED=
SKIP_MAPPING=false
SKIP_POOLING=false
SKIP_CALLING=false
while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -b | --barcode-length ) BARCODE_LEN="$2"; shift 2 ;;
    -e | --exp-prefix ) EXP_PREFIX="$2"; shift 2 ;;
    -f | --fasta ) FASTA="$2"; shift 2 ;;
    -g | --gsize ) GENOME_SIZE="$2"; shift 2 ;;
    -i | --input ) INPUT_DIR="$2"; shift 2 ;;
    -m | --min-length ) MIN_LENGTH="$2"; shift 2 ;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2 ;;
    -p | --pool-prefix ) POOL_PREFIX="$2"; shift 2 ;;
    -q | --quality ) QUALITY="$2"; shift 2 ;;
    -t | --thread ) THREAD="$2"; shift 2 ;;
    --dbkey ) DB_KEY="$2"; shift 2 ;;
    --full-bed ) FULL_BED="$2"; shift 2 ;;
    --index ) BWA_INDEX="$2"; shift 2 ;;
    --longest-bed ) LONGEST_BED="$2"; shift 2 ;;
    --repeat-bed ) REPEAT_BED="$2"; shift 2 ;;
    --skip-mapping ) SKIP_MAPPING=true; shift ;;
    --skip-pooling ) SKIP_POOLING=true; shift ;;
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

if [ ! -f ${BWA_INDEX}.bwt ]; then
  echo "--index: Wrong! Invalid BWA index!"
  exit 2
fi

# basic variables
REP_NUM=`find $INPUT_DIR -type f -name "${EXP_PREFIX}*.trim.fastq" | wc -l`
MAP_THREAD=$((THREAD / REP_NUM))

# basic variables
FILT_DIR="$OUTPUT_DIR/filter"
MAP_DIR="$OUTPUT_DIR/mapping"
FASTQC_DIR="$OUTPUT_DIR/mapping/fastQC"
CLUSTER_DIR="$OUTPUT_DIR/cluster"
CIMS_DIR="$OUTPUT_DIR/CIMS"
CITS_DIR="$OUTPUT_DIR/CITS"
FINAL_DIR="$OUTPUT_DIR/final"
POOL_PREFIX="HepG2_m6A_WT_IP"
MAX_DIFF=3

if [[ ! -d $OUTPUT_DIR ]]; then
  mkdir -p $OUTPUT_DIR
fi

if [[ ! -d $FILT_DIR  ]]; then
    mkdir -p $FILT_DIR
fi

if [[ ! -d $MAP_DIR  ]]; then
    mkdir -p $MAP_DIR
fi

if $SKIP_MAPPING; then
  echo "Skip mapping."
else
  # Estimating reads length
  echo "Estimating reads length..."
  READ_AVG_LENGTH=100
  for i in "$INPUT_DIR/${EXP_PREFIX}*.fastq"; do
    READ_AVG_LENGTH=`awk '{if(NR%4==2) {count++; lengthSum += length} } END{print int(lengthSum/count)}' $i`
    break 1
  done
  echo "Reads length:${READ_AVG_LENGTH}..."
  
  if (( $READ_AVG_LENGTH <= 17 )); then
    MAX_DIFF=1
  elif (( $READ_AVG_LENGTH <= 20 )); then
    MAX_DIFF=2
  elif (( $READ_AVG_LENGTH <= 45 )); then
    MAX_DIFF=3
  elif (( $READ_AVG_LENGTH <= 73 )); then
    MAX_DIFF=4
  elif (( $READ_AVG_LENGTH <= 104 )); then
    MAX_DIFF=5
  elif (( $READ_AVG_LENGTH <= 137 )); then
    MAX_DIFF=6
  elif (( $READ_AVG_LENGTH <= 172 )); then
    MAX_DIFF=7
  elif (( $READ_AVG_LENGTH <= 208 )); then
    MAX_DIFF=8
  elif (( $READ_AVG_LENGTH <= 244 )); then
    MAX_DIFF=9
  fi
  echo "max_diff in BWA:${MAX_DIFF}"
  # Read preprocessing
  ## parallel start
  echo "Mapping with BWA."
  REP_NUM=`find $INPUT_DIR -type f -name "${EXP_PREFIX}*.trim.fastq" | wc -l`
  MAP_THREAD=$((THREAD / REP_NUM))
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
      cd $FILT_DIR
      PREFIX=${i%%.fastq}
      PREFIX=${PREFIX##*/} # *.trim
      MAP_PREFIX=${i%%.trim.fastq}
      MAP_PREFIX=${MAP_PREFIX##*/}
      fastq2collapse.pl $i - | gzip -c > ${PREFIX}.c.fastq.gz
      
      ## Strip random barcode (UMI)
      stripBarcode.pl -format fastq -len ${BARCODE_LEN} ${PREFIX}.c.fastq.gz - | \
        gzip -c > ${PREFIX}.c.tag.fastq.gz
      
      zcat -f ${PREFIX}.c.tag.fastq.gz | awk '{if(NR%4==2) {print length($0)}}' | \
        sort -n | uniq -c | awk '{print $2"\t"$1}' \
        > ${PREFIX}.c.tag.seqlen.stat.txt
      
      ## Read mapping & parsing
      cd $MAP_DIR
      
      ## Read $MAP_DIR
      bwa aln -t ${MAP_THREAD} -n ${MAX_DIFF} -q ${QUALITY} \
        ${BWA_INDEX} $FILT_DIR/${PREFIX}.c.tag.fastq.gz \
        > ${PREFIX}.sai 2> ${PREFIX}.bwa.log
      bwa samse $BWA_INDEX ${PREFIX}.sai $FILT_DIR/${PREFIX}.c.tag.fastq.gz | \
        gzip -c > ${MAP_PREFIX}.sam.gz
      
      ## Parsing SAM file
      parseAlignment.pl -v --map-qual 1 --min-len $MIN_LENGTH --mutation-file \
        ${MAP_PREFIX}.mutation.txt ${MAP_PREFIX}.sam.gz ${MAP_PREFIX}.tag.bed \
        > ${MAP_PREFIX}.parseAlignment.log 2>&1
  
      # Remove tags from rRNA and other repetitive RNA (optional)
      tagoverlap.pl -big -region ${REPEAT_BED} -ss --complete-overlap \
        -r --keep-tag-name --keep-score -v ${MAP_PREFIX}.tag.bed \
        ${MAP_PREFIX}.tag.norRNA.bed > ${MAP_PREFIX}.tagoverlap.log 2>&1
      # Collapse PCR duplicates
      tag2collapse.pl -big -v --random-barcode -EM 30 \
        --seq-error-model alignment -weight --weight-in-name \
        --keep-max-score --keep-tag-name ${MAP_PREFIX}.tag.norRNA.bed \
        ${MAP_PREFIX}.tag.uniq.bed
      awk '{print $3-$2}' ${MAP_PREFIX}.tag.uniq.bed | sort -n | uniq -c | \
        awk '{print $2"\t"$1}' > ${MAP_PREFIX}.tag.uniq.len.dist.txt
      
      ## Get the mutations in unique tags
      joinWrapper.py ${MAP_PREFIX}.mutation.txt \
        ${MAP_PREFIX}.tag.uniq.bed 4 4 N ${MAP_PREFIX}.tag.uniq.mutation.txt
      ## end
      echo >&9
    } &
  done
  
  wait
  exec 9>&-
  exec 9<&-
  ## parallel end
fi

## Merging biological replicates
if $SKIP_POOLING; then
  echo "Skip pooling."
else
  cd $MAP_DIR
  declare -a colorArr
  colorArr=("153,0,0" "153,76,0" "76,153,0" "0,76,153" "76,0,153" "153,76,0" "0,153,153")
  for ((i=1; i<=$REP_NUM; i++))
  do
    rep=$i
    colorIndex=$((i-1))
    color="${colorArr[$colorIndex]}"
    bed2rgb.pl -v -col $color ${EXP_PREFIX}${rep}.tag.uniq.bed ${EXP_PREFIX}${rep}.tag.uniq.rgb.bed
  done
  
  cat ${EXP_PREFIX}*.tag.uniq.rgb.bed > ${POOL_PREFIX}.pool.tag.uniq.rgb.bed
  cat ${EXP_PREFIX}*.tag.uniq.mutation.txt > ${POOL_PREFIX}.pool.tag.uniq.mutation.txt
  
  ## Annotating and visualizing CLIP tags
  bed2annotation.pl -dbkey ${DB_KEY} -ss -big -region -v \
    -summary ${POOL_PREFIX}.pool.tag.uniq.annot.summary.txt \
    ${POOL_PREFIX}.pool.tag.uniq.rgb.bed \
    ${POOL_PREFIX}.pool.tag.uniq.annot.txt \
    > ${POOL_PREFIX}.tag.annot.log 2>&1
  
  ## Generate bedgraph for visualization in the genome browser
  tag2profile.pl -v -ss -exact -of bedgraph -n "${POOL_PREFIX} Unique Tag Profile" \
    ${POOL_PREFIX}.pool.tag.uniq.rgb.bed \
    ${POOL_PREFIX}.pool.tag.uniq.bedgraph \
    > ${POOL_PREFIX}.tag.bedgraph.log 2>&1
fi


# Peak calling, cluster
if [[ ! -d $CLUSTER_DIR  ]]; then
    mkdir $CLUSTER_DIR
fi

cd $CLUSTER_DIR

ln -sf $MAP_DIR/${POOL_PREFIX}.pool.tag.uniq.rgb.bed $CLUSTER_DIR

if $SKIP_CALLING; then
  echo "Skip peak calling."
else
  ## Mode 1: Peak calling with no statistical significance
  PEAK_PREFIX="${POOL_PREFIX}.peak.nostats"
  tag2peak.pl -big -ss -v --valley-seeking --valley-depth 0.9 \
    ${POOL_PREFIX}.pool.tag.uniq.rgb.bed ${PEAK_PREFIX}.bed \
    --out-boundary ${PEAK_PREFIX}.boundary.bed \
    --out-half-PH ${PEAK_PREFIX}.halfPH.bed \
    > ${PEAK_PREFIX}.mode1.log 2>&1
  
  bed2annotation.pl -dbkey ${DB_KEY} -ss -big -region -v \
    -summary ${PEAK_PREFIX}.annot.summary.txt \
    ${PEAK_PREFIX}.bed \
    ${PEAK_PREFIX}.annot.txt \
    > ${PEAK_PREFIX}.annot.log 2>&1
  ### scan RRACH motif
  scanMotif.py -input ${PEAK_PREFIX}.bed -format bed6 \
    -fasta ${FASTA} -motif RRACH -tag 3 \
    -output ${PEAK_PREFIX}.RRACH.bed
  
  ## Mode 2: Peak calling with statistical significance
  PEAK_PREFIX="${POOL_PREFIX}.peak.sig"
  tag2peak.pl -big -ss -v --valley-seeking -p 0.05 \
    --valley-depth 0.9 --multi-test --dbkey ${DB_KEY} \
    ${POOL_PREFIX}.pool.tag.uniq.rgb.bed ${PEAK_PREFIX}.bed \
    --out-boundary ${PEAK_PREFIX}.boundary.bed \
    --out-half-PH ${PEAK_PREFIX}.halfPH.bed \
    > ${PEAK_PREFIX}.mode2.log 2>&1
  
  ### scan RRACH motif
  scanMotif.py -input ${PEAK_PREFIX}.bed -format bed6 \
    -fasta ${FASTA} -motif RRACH -tag 3 \
    -output ${PEAK_PREFIX}.RRACH.bed
fi

# CIMS calling
if [[ ! -d $CIMS_DIR  ]]; then
    mkdir $CIMS_DIR
fi
cd $CIMS_DIR
ln -sf $MAP_DIR/${POOL_PREFIX}.pool.tag.uniq.rgb.bed $CIMS_DIR
ln -sf $MAP_DIR/${POOL_PREFIX}.pool.tag.uniq.mutation.txt $CIMS_DIR

if $SKIP_CALLING; then
  echo "Skip CIMS calling."
else
  ## Get specific types of mutations
  ### del: deletions, ins: insertions, sub: substitutions
  echo "CIMS calling..."
  getMutationType.pl -t del ${POOL_PREFIX}.pool.tag.uniq.mutation.txt ${POOL_PREFIX}.del.bed
  getMutationType.pl -t ins ${POOL_PREFIX}.pool.tag.uniq.mutation.txt ${POOL_PREFIX}.ins.bed
  getMutationType.pl -t sub ${POOL_PREFIX}.pool.tag.uniq.mutation.txt ${POOL_PREFIX}.sub.bed
  
  CIMS.pl -big -n 10 -p -outp ${POOL_PREFIX}.sub.posStat.txt \
    -v ${POOL_PREFIX}.pool.tag.uniq.rgb.bed \
    ${POOL_PREFIX}.sub.bed \
    ${POOL_PREFIX}.sub.CIMS.txt \
    > ${POOL_PREFIX}.sub.CIMS.log 2>&1
fi
sort -k 9,9n -k 8,8nr -k 7,7n ${POOL_PREFIX}.sub.CIMS.txt | \
  cut -f 1-6 > ${POOL_PREFIX}.sub.CIMS.bed

## get C->T mutations
ctk_C2T_mutation_filter.sh "${POOL_PREFIX}.sub.CIMS.bed" \
  "${POOL_PREFIX}.pool.tag.uniq.mutation.txt" \
  "${POOL_PREFIX}.sub.CIMS.CT.bed"

## get m6A sites with RRACH motif
bedtools shift -i ${POOL_PREFIX}.sub.CIMS.CT.bed \
  -g ${GENOME_SIZE} -m 1 -p -1 | \
  bedtools slop -i stdin -b 2 -g ${GENOME_SIZE} \
  > ${POOL_PREFIX}.temp.bed

scanMotif.py -input ${POOL_PREFIX}.temp.bed -format bed6 \
  -fasta ${FASTA} -motif RRACH -tag 3 \
  -output ${POOL_PREFIX}.sub.CIMS.CT.RRACH.bed

rm ${POOL_PREFIX}.temp.bed

##significant: get CIMS
awk '{if($9<=0.05) {print $0}}' ${POOL_PREFIX}.sub.CIMS.txt | \
  sort -k 9,9n -k 8,8nr -k 7,7n > ${POOL_PREFIX}.sub.CIMS.sig.txt
cut -f 1-6 ${POOL_PREFIX}.sub.CIMS.sig.txt > ${POOL_PREFIX}.sub.CIMS.sig.bed

##significant: get C->T mutations
ctk_C2T_mutation_filter.sh "${POOL_PREFIX}.sub.CIMS.sig.bed" \
  "${POOL_PREFIX}.pool.tag.uniq.mutation.txt" \
  "${POOL_PREFIX}.sub.CIMS.sig.CT.bed"

##significant: get m6A sites with RRACH motif
bedtools shift -i ${POOL_PREFIX}.sub.CIMS.sig.CT.bed \
  -g ${GENOME_SIZE} -m 1 -p -1 | \
  bedtools slop -i stdin -b 2 -g ${GENOME_SIZE} \
  > ${POOL_PREFIX}.temp.bed

scanMotif.py -input ${POOL_PREFIX}.temp.bed -format bed6 \
  -fasta ${FASTA} -motif RRACH -tag 3 \
  -output ${POOL_PREFIX}.sub.CIMS.sig.CT.RRACH.bed

rm ${POOL_PREFIX}.temp.bed

# CITS calling
if [[ ! -d $CITS_DIR  ]]; then
    mkdir $CITS_DIR
fi
cd $CITS_DIR

ln -sf $MAP_DIR/${POOL_PREFIX}.pool.tag.uniq.rgb.bed
ln -sf $CIMS_DIR/${POOL_PREFIX}.del.bed

if $SKIP_CALLING; then
  echo "Skip CITS calling."
else
  ## all CITS
  echo "CITS calling..."
  CITS.pl -big -p 0.999 --gap 25 -v \
    ${POOL_PREFIX}.pool.tag.uniq.rgb.bed \
    ${POOL_PREFIX}.del.bed \
    ${POOL_PREFIX}.CITS.bed \
    > ${POOL_PREFIX}.CITS.log 2>&1
  ## significant CITS
  CITS.pl -big -p 0.05 --gap 25 -v \
    ${POOL_PREFIX}.pool.tag.uniq.rgb.bed \
    ${POOL_PREFIX}.del.bed \
    ${POOL_PREFIX}.CITS.sig.bed \
    > ${POOL_PREFIX}.CITS.sig.log 2>&1
fi

bedtools shift -i ${POOL_PREFIX}.CITS.bed \
  -g ${GENOME_SIZE} -m 1 -p -1 | \
  bedtools slop -i stdin -b 2 -g ${GENOME_SIZE} \
  > ${POOL_PREFIX}.temp.bed

scanMotif.py -input ${POOL_PREFIX}.temp.bed -format bed6 \
  -fasta ${FASTA} -motif RRACH -tag 3 \
  -output ${POOL_PREFIX}.CITS.RRACH.bed

## get m6A sites with RRACH motif
bedtools shift -i ${POOL_PREFIX}.CITS.sig.bed \
  -g ${GENOME_SIZE} -m 1 -p -1 | \
  bedtools slop -i stdin -b 2 -g ${GENOME_SIZE} \
  > ${POOL_PREFIX}.temp.bed

scanMotif.py -input ${POOL_PREFIX}.temp.bed -format bed6 \
  -fasta ${FASTA} -motif RRACH -tag 3 \
  -output ${POOL_PREFIX}.CITS.sig.RRACH.bed

rm ${POOL_PREFIX}.temp.bed

# generating the final results
if [[ ! -d $FINAL_DIR  ]]; then
    mkdir $FINAL_DIR
fi
cd $FINAL_DIR

cp $CIMS_DIR/${POOL_PREFIX}.sub.CIMS.CT.RRACH.bed ./
cp $CIMS_DIR/${POOL_PREFIX}.sub.CIMS.sig.CT.RRACH.bed ./
cp $CITS_DIR/${POOL_PREFIX}.CITS.RRACH.bed ./
cp $CITS_DIR/${POOL_PREFIX}.CITS.sig.RRACH.bed ./

cat ${POOL_PREFIX}.sub.CIMS.CT.RRACH.bed ${POOL_PREFIX}.CITS.RRACH.bed | \
  sort -t $'\t' -k1,1 -k2,2n > ${POOL_PREFIX}.combine.RRACH.bed
cat ${POOL_PREFIX}.sub.CIMS.sig.CT.RRACH.bed ${POOL_PREFIX}.CITS.sig.RRACH.bed | \
  sort -t $'\t' -k1,1 -k2,2n > ${POOL_PREFIX}.combine.sig.RRACH.bed

if [ ! -z $LONGEST_BED ]; then
  for i in `find ./ -type f -name "${POOL_PREFIX}*.bed"`;
  do
    PREFIX=${i%%.bed}
    bedBinDistribution.pl -input $i -bed12 $LONGEST_BED \
      -o ${PREFIX}.percentage.bin
    bedBinDistribution.pl -input $i -bed12 $LONGEST_BED \
      --type count -o ${PREFIX}.count.bin
  done
fi

if [ ! -z $FULL_BED ]; then
  ## get mRNA annotation bed12
  cat $FULL_BED | awk '/protein_coding.+protein_coding\t/' > mRNA.annotation.bed12.tmp
  for i in `find ./ -type f -name "${POOL_PREFIX}*.bed"`;
  do
    PREFIX=${i%%.bed}
    ### gene type
    geneDistribution.pl -strand --input $i \
      -bed12 $FULL_BED -o ${PREFIX}.gene
    sed -i '1i geneType\tpeakNumber' ${PREFIX}.gene
    ### gene region
    regionDistribution.pl -strand -size 200 -f '5utr,cds,stopCodon,3utr' \
      --input $i \
      -bed12 mRNA.annotation.bed12.tmp -o ${PREFIX}.region
    sed -i '1i region\tpeakNumber\tenrichment' ${PREFIX}.region
  done
  rm -f mRNA.annotation.bed12.tmp
fi

rm -f *.tmp
