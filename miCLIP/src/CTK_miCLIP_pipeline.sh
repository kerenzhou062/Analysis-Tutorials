#!/bin/bash
#SBATCH --job-name=CTK_miCLIP_pipeline    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail  
#SBATCH -n 10                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=100G                      # Amount of memory in GB
#SBATCH --time=120:10:00               # Time limit hrs:min:sec
#SBATCH --output=CTK_miCLIP_pipeline.log   # Standard output and error log

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.

function showHelp {
  echo -ne "usage: sbatch CTK_miCLIP_pipeline.sh -n <thread_num> -o <log> --mem <200G> "
  echo -ne "CTK_miCLIP_pipeline.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
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
    --mfreq: # of mutation tags (defined enriched sites) <int>
    --mkr: m/k ratio for enriched mutation sites <str>
    --dbkey: CTK dbkey (hg38|hg19|mm10) <str>
    --gene: custom gene annotation bed (tag2peak.pl) <str>
    --full-bed: all transcripts annotation in bed12 <str>
    --index: genome index <str>
    --joinWrapper-loc: # location of joinWrapper.py <str>
    --longest-bed: mRNA longest annotation in bed12 <str>
    --motif: MOTIF sequence used to search (RRACH) <str>
    --mtag: -tag parameter in scanMotif.py <str>
    --repeat-bed: repeat bed used for filtering (eg.t/rRNA) <str>
    --keep-tmp-fastq: keep temporary fastqs <bool>
    --skip-mapping: skip reads mapping step <bool>
    --skip-pooling: skip pooling step <bool>
    --skip-PCR-am: skip removing PCR duplicates after mapping <bool>
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
  --long keep-tmp-fastq,skip-PCR-am \
  --long input:,thread:,output:,exp-prefix:,pool-prefix:,index:,min-length:,dbkey:,joinWrapper-loc: \
  --long longest-bed:,full-bed:,repeat-bed:,barcode-length:,fasta:,quality:,mkr:,motif:,mtag: \
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
CUSTOM_GENE=
QUALITY=20
FASTA=
GENOME_SIZE=
LONGEST_BED=
JOIN_WRAPPER_LOC=
FULL_BED=
REPEAT_BED=
MUTATE_FREQ=1
MKR_RATIO=0.5
MOTIF="RRACH"
MOTIF_TAG=3
KEEP_TMP_FASTQ=false
SKIP_MAPPING=false
SKIP_POOLING=false
SKIP_PCR_AM=false
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
    --mfreq ) MUTATE_FREQ="$2"; shift 2 ;;
    --mkr ) MKR_RATIO="$2"; shift 2 ;;
    --dbkey ) DB_KEY="$2"; shift 2 ;;
    --gene ) CUSTOM_GENE="$2"; shift 2 ;;
    --full-bed ) FULL_BED="$2"; shift 2 ;;
    --index ) BWA_INDEX="$2"; shift 2 ;;
    --joinWrapper-loc ) JOIN_WRAPPER_LOC="$2"; shift 2 ;;
    --longest-bed ) LONGEST_BED="$2"; shift 2 ;;
    --motif ) MOTIF="$2"; shift 2 ;;
    --mtag ) MOTIF_TAG="$2"; shift 2 ;;
    --repeat-bed ) REPEAT_BED="$2"; shift 2 ;;
    --keep-tmp-fastq ) KEEP_TMP_FASTQ=true; shift ;;
    --skip-mapping ) SKIP_MAPPING=true; shift ;;
    --skip-pooling ) SKIP_POOLING=true; shift ;;
    --skip-PCR-am ) SKIP_PCR_AM=true; shift ;;
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

if [ -z $DB_KEY ]; then
  echo "--dbkey: Wrong! Please specify --dbkey";
fi

if [ ! -f ${BWA_INDEX}.bwt ]; then
  echo "--index: Wrong! Invalid BWA index!"
  exit 2
fi

echo "Running pipleline with following parameters:"
echo "THREAD=$THREAD"
echo "BARCODE_LEN=$BARCODE_LEN"
echo "MIN_LENGTH=$MIN_LENGTH"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo "INPUT_DIR=$INPUT_DIR"
echo "EXP_PREFIX=$EXP_PREFIX"
echo "POOL_PREFIX=$POOL_PREFIX"
echo "BWA_INDEX=$BWA_INDEX"
echo "DB_KEY=$DB_KEY"
echo "CUSTOM_GENE=$CUSTOM_GENE"
echo "QUALITY=$QUALITY"
echo "FASTA=$FASTA"
echo "GENOME_SIZE=$GENOME_SIZE"
echo "LONGEST_BED=$LONGEST_BED"
echo "JOIN_WRAPPER_LOC=$JOIN_WRAPPER_LOC"
echo "FULL_BED=$FULL_BED"
echo "REPEAT_BED=$REPEAT_BED"
echo "MUTATE_FREQ=$MUTATE_FREQ"
echo "MKR_RATIO=$MKR_RATIO"
echo "MOTIF=$MOTIF"
echo "MOTIF_TAG=$MOTIF_TAG"
echo "KEEP_TMP_FASTQ=$KEEP_TMP_FASTQ"
echo "SKIP_MAPPING=$SKIP_MAPPING"
echo "SKIP_POOLING=$SKIP_POOLING"
echo "SKIP_PCR_AM=$SKIP_PCR_AM"
echo "SKIP_CALLING=$SKIP_CALLING"
echo ""

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
MAX_DIFF=3
SLOP_B=$((${#MOTIF}/2))

if [[ ! -d $OUTPUT_DIR ]]; then
  mkdir -p $OUTPUT_DIR
fi

if [[ ! -d $FILT_DIR  ]]; then
  mkdir -p $FILT_DIR
fi

if [[ ! -d $MAP_DIR  ]]; then
  mkdir -p $MAP_DIR
fi

if [[ -z $JOIN_WRAPPER_LOC ]]; then
  joinWrapper="joinWrapper.py"
else
  joinWrapper="$JOIN_WRAPPER_LOC/joinWrapper.py"
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
    {
      cd $FILT_DIR
      PREFIX=${i%%.fastq}
      PREFIX=${PREFIX##*/} # *.trim
      MAP_PREFIX=${i%%.trim.fastq}
      MAP_PREFIX=${MAP_PREFIX##*/}

      ## Collapse exact duplicates
      fastq2collapse.pl $i - | gzip -c > ${PREFIX}.c.fastq.gz
      
      if (( $BARCODE_LEN > 0 )); then
        ## Strip random barcode (UMI)
        stripBarcode.pl -format fastq -len ${BARCODE_LEN} ${PREFIX}.c.fastq.gz - | \
          gzip -c > ${PREFIX}.c.tag.fastq.gz
        zcat -f ${PREFIX}.c.tag.fastq.gz | awk '{if(NR%4==2) {print length($0)}}' | \
          sort -n | uniq -c | awk '{print $2"\t"$1}' \
          > ${PREFIX}.c.tag.seqlen.stat.txt
        input=${PREFIX}.c.tag.fastq.gz
      else
        input=${PREFIX}.c.fastq.gz
      fi
      
      ## Read mapping & parsing
      cd $MAP_DIR
      
      ## Read $MAP_DIR
      bwa aln -t ${MAP_THREAD} -n ${MAX_DIFF} -q ${QUALITY} \
        ${BWA_INDEX} $FILT_DIR/${input} \
        > ${PREFIX}.sai 2> ${PREFIX}.bwa.log
      bwa samse $BWA_INDEX ${PREFIX}.sai $FILT_DIR/${input} | \
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
      if (( $BARCODE_LEN > 0 )); then
        if $SKIP_PCR_AM; then
          cp ${MAP_PREFIX}.tag.norRNA.bed ${MAP_PREFIX}.tag.uniq.bed
        else
          tag2collapse.pl -big -v --random-barcode -EM 30 \
            --seq-error-model alignment -weight --weight-in-name \
            --keep-max-score --keep-tag-name ${MAP_PREFIX}.tag.norRNA.bed \
            ${MAP_PREFIX}.tag.uniq.bed
        fi
      else
        if $SKIP_PCR_AM; then
          cp ${MAP_PREFIX}.tag.norRNA.bed ${MAP_PREFIX}.tag.uniq.bed
        else
          tag2collapse.pl -v -big -weight --weight-in-name --keep-max-score \
            --keep-tag-name ${MAP_PREFIX}.tag.norRNA.bed ${MAP_PREFIX}.tag.uniq.bed
        fi
      fi
      
      awk '{print $3-$2}' ${MAP_PREFIX}.tag.uniq.bed | sort -n | uniq -c | \
        awk '{print $2"\t"$1}' > ${MAP_PREFIX}.tag.uniq.len.dist.txt
      
      ## Get the mutations in unique tags
      $joinWrapper ${MAP_PREFIX}.mutation.txt \
        ${MAP_PREFIX}.tag.uniq.bed 4 4 N ${MAP_PREFIX}.tag.uniq.mutation.txt
      ## end
      echo >&9
    } &
  done
  
  wait
  exec 9>&-
  exec 9<&-
  ## parallel end
  #delete tmp fastqs
  if $KEEP_TMP_FASTQ; then
    echo "Keep tmp fastqs."
  else
    echo "Deleting tmp fastqs..."
    find $FILT_DIR -maxdepth 1 -type f -name "${EXP_PREFIX}*.fastq" | xargs -I {} rm -f {}
    find $FILT_DIR -maxdepth 1 -type f -name "${EXP_PREFIX}*.fastq.gz" | xargs -I {} rm -f {}
  fi
fi

## Merging biological replicates
if $SKIP_POOLING; then
  echo "Skip pooling."
else
  echo $REP;
  cd $MAP_DIR
  declare -a colorArr
  colorArr=("153,0,0" "153,76,0" "76,153,0" "0,76,153" "76,0,153" "153,76,0" "0,153,153")
  for ((i=1; i<=$REP_NUM; i++))
  do
    rep=$i
    colorIndex=$((i-1))
    color="${colorArr[$colorIndex]}"
    if (( $REP > 1 )); then
      bed2rgb.pl -v -col $color ${EXP_PREFIX}${rep}.tag.uniq.bed ${EXP_PREFIX}${rep}.tag.uniq.rgb.bed
    else
      bed2rgb.pl -v -col $color ${EXP_PREFIX}.tag.uniq.bed ${EXP_PREFIX}.tag.uniq.rgb.bed
    fi
  done

  if (( $REP > 1 )); then
    cp ${EXP_PREFIX}.tag.uniq.rgb.bed ${POOL_PREFIX}.pool.tag.uniq.rgb.bed
  else
    cat ${EXP_PREFIX}*.tag.uniq.rgb.bed > ${POOL_PREFIX}.pool.tag.uniq.rgb.bed
    cat ${EXP_PREFIX}*.tag.uniq.mutation.txt > ${POOL_PREFIX}.pool.tag.uniq.mutation.txt
  fi
  
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
  echo "Peak calling [mode1]..."
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
    -fasta ${FASTA} -motif RRACH -tag ${MOTIF_TAG} \
    -output ${PEAK_PREFIX}.${MOTIF}.bed
  
  ## Mode 2: Peak calling with statistical significance
  echo "Peak calling [mode2]..."
  PEAK_PREFIX="${POOL_PREFIX}.peak.sig"
  if [[ -z $CUSTOM_GENE ]]; then
    tag2peak.pl -big -ss -v --valley-seeking -p 0.05 \
      --valley-depth 0.9 --multi-test --dbkey ${DB_KEY} \
      ${POOL_PREFIX}.pool.tag.uniq.rgb.bed ${PEAK_PREFIX}.bed \
      --out-boundary ${PEAK_PREFIX}.boundary.bed \
      --out-half-PH ${PEAK_PREFIX}.halfPH.bed \
      > ${PEAK_PREFIX}.mode2.log 2>&1
  else
    echo "Peak calling with CUSTOM GENE annotation..."
    tag2peak.pl -big -ss -v --valley-seeking -p 0.05 \
      --valley-depth 0.9 --multi-test --gene ${CUSTOM_GENE} \
      ${POOL_PREFIX}.pool.tag.uniq.rgb.bed ${PEAK_PREFIX}.bed \
      --out-boundary ${PEAK_PREFIX}.boundary.bed \
      --out-half-PH ${PEAK_PREFIX}.halfPH.bed \
      > ${PEAK_PREFIX}.mode2.log 2>&1
  fi
  
  ### scan RRACH motif
  scanMotif.py -input ${PEAK_PREFIX}.bed -format bed6 \
    -fasta ${FASTA} -motif RRACH -tag ${MOTIF_TAG} \
    -output ${PEAK_PREFIX}.${MOTIF}.bed
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
  getMutationType.pl -t c2t ${POOL_PREFIX}.pool.tag.uniq.mutation.txt ${POOL_PREFIX}.c2t.bed

  CIMS.pl -big -n 10 -p -outp ${POOL_PREFIX}.c2t.posStat.txt \
    -v ${POOL_PREFIX}.pool.tag.uniq.rgb.bed \
    ${POOL_PREFIX}.c2t.bed \
    ${POOL_PREFIX}.c2t.CIMS.txt \
    > ${POOL_PREFIX}.c2t.CIMS.log 2>&1
fi

## calculating m/k ratio: mutationFreq(m)/tagNumber(k)
awk 'BEGIN{FS="\t";OFS="\t";}{if(FNR==1){$5=$5"(m/k)"; print $0}else{$5=$8/$7;print $0}}' \
  ${POOL_PREFIX}.c2t.CIMS.txt > ${POOL_PREFIX}.c2t.CIMS.mk.txt

awk 'BEGIN{FS="\t";OFS="\t";}{if(FNR>1){print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' \
  ${POOL_PREFIX}.c2t.CIMS.mk.txt | sort -k 5,5nr -k 8,8nr -k 7,7n \
  > ${POOL_PREFIX}.c2t.CIMS.bed

##enriched: get CIMS, m/k>=0.5
awk -v mkr="${MKR_RATIO}" -v mfreq=${MUTATE_FREQ} \
  'BEGIN{FS="\t";OFS="\t";}{if($5>=mkr && $8>=mfreq) {print $0}}' \
  ${POOL_PREFIX}.c2t.CIMS.mk.txt > ${POOL_PREFIX}.c2t.CIMS.enrich.bed

##significant: get CIMS
awk 'BEGIN{FS="\t";OFS="\t";}{if(FNR==1){print $0}else{if($9<=0.05) {print $0}}}' \
  ${POOL_PREFIX}.c2t.CIMS.mk.txt > ${POOL_PREFIX}.c2t.CIMS.sig.txt

awk 'BEGIN{FS="\t";OFS="\t";}{if(FNR>1){print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' \
  ${POOL_PREFIX}.c2t.CIMS.sig.txt | sort -k 5,5nr -k 8,8nr -k 7,7n \
  > ${POOL_PREFIX}.c2t.CIMS.sig.bed

##get m6A sites with ${MOTIF} motif
for i in `find ./ -type f -name "${POOL_PREFIX}*.bed" | grep "CIMS" | grep -v "${MOTIF}"`;
do
  prefix=${i%%.bed}
  sort -t $'\t' -k1,1n -k2,2n ${i} | bedtools shift -i stdin -g ${GENOME_SIZE} -p -1 -m 1 | \
    bedtools slop -i stdin -b $SLOP_B -s -g ${GENOME_SIZE} > ${POOL_PREFIX}.temp.bed
  scanMotif.py -input ${POOL_PREFIX}.temp.bed -format bed6 \
    -fasta ${FASTA} -motif ${MOTIF} -tag ${MOTIF_TAG} \
    -output ${prefix}.${MOTIF}.bed
done

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

##get m6A sites with ${MOTIF} motif
for i in `find ./ -type f -name "${POOL_PREFIX}*.bed" | grep "CITS" | grep -v "${MOTIF}"`;
do
  prefix=${i%%.bed}
  sort -t $'\t' -k1,1n -k2,2n ${i} | bedtools shift -i stdin -g ${GENOME_SIZE} -p -1 -m 1 | \
    bedtools slop -i stdin -b $SLOP_B -s -g ${GENOME_SIZE} > ${POOL_PREFIX}.temp.bed
  scanMotif.py -input ${POOL_PREFIX}.temp.bed -format bed6 \
    -fasta ${FASTA} -motif ${MOTIF} -tag ${MOTIF_TAG} \
    -output ${prefix}.${MOTIF}.bed
done

rm ${POOL_PREFIX}.temp.bed

# generating the final results
if [[ ! -d $FINAL_DIR  ]]; then
  mkdir $FINAL_DIR
fi
cd $FINAL_DIR

cp $CIMS_DIR/${POOL_PREFIX}.c2t.CIMS.${MOTIF}.bed ./
cp $CIMS_DIR/${POOL_PREFIX}.c2t.CIMS.enrich.${MOTIF}.bed ./
cp $CIMS_DIR/${POOL_PREFIX}.c2t.CIMS.sig.${MOTIF}.bed ./
cp $CITS_DIR/${POOL_PREFIX}.CITS.${MOTIF}.bed ./
cp $CITS_DIR/${POOL_PREFIX}.CITS.sig.${MOTIF}.bed ./

## cat CIMS and CITS and uniq
function combine {
  input1=$1
  input2=$2
  output=$3
  cat $input1 $input2 | \
    awk 'BEGIN{OFS="\t";FS="\t";}
    {
      if($1!="#chrom") {
        key = $1"\t"$2"\t"$3":"$5"\t"$6;
        if (match($4, /CITS/) ) {
          type = "CITS";
        }else{
          type = "CIMS";
        }
        if (key in arrayA) {
          arrayA[key] = arrayA[key]"|"type
        }else{
          arrayA[key] = $(NF-1)"="type
        }
      }
    }
    END{
      for (key in arrayA) {
        split(key,splitArr,":");
        pos = splitArr[1];
        remain = splitArr[2];
        print pos, arrayA[key], remain;
      }
    }' | sort -k1,1 -k2,2n > $output
}

combine ${POOL_PREFIX}.c2t.CIMS.${MOTIF}.bed \
  ${POOL_PREFIX}.CITS.${MOTIF}.bed \
  ${POOL_PREFIX}.combine.${MOTIF}.bed
combine ${POOL_PREFIX}.c2t.CIMS.sig.${MOTIF}.bed \
  ${POOL_PREFIX}.CITS.sig.${MOTIF}.bed \
  ${POOL_PREFIX}.combine.sig.${MOTIF}.bed

echo "Pooling CT and Truncation m6A sites done."

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
rm -f *.tmp

echo "beds annotation done."
