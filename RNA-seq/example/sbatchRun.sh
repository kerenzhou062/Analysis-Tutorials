#/bin/sh
PIPELINE_DIR=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/home/public/pipeline/RNA-seq

BASE=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/publicData/
LOG_DIR="$BASE/log"
FASTQ_DIR="$BASE/fastq"
ALIGN_DIR="$BASE/alignment"
SAMPLE_MATRIX="$BASE/sample.matrix"
# rsem reference
PUBLIC_BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome/index/rsem/encode"
STAR_GENOME_DIR="$PUBLIC_BASE/star/hg38/v32_tsnoRNA"
RSEM_GENOME_DIR="$PUBLIC_BASE/rsem/hg38/v32_tsnoRNA/RSEMref"
CPU=8
MEMORY=80


cd $BASE

function runRSEM {
  EXP_PREFIX=$1
  READ1=$2
  READ2=$3
  SEQ_TYPE=$4
  STRAND=$5
  TEST=$6
  if [[ $READ2 == "none" ]]; then
    READ2=""
  fi
  if [[ $TEST == "test" ]]; then
    # run RSEM_STAR_align_pipeline
    sbatch -n $CPU --mem ${MEMORY}G -o $LOG_DIR/RSEM_STAR_align_pipeline.STAR.${EXP_PREFIX}.log -J RSEM_STAR_align_pipeline_${EXP_PREFIX} \
      $PIPELINE_DIR/RSEM_STAR_align_pipeline.sh --thread $CPU --read1 $READ1 --read2 $READ2 --mem $MEMORY -o $ALIGN_DIR/$EXP_PREFIX \
      -p ${EXP_PREFIX} --seq-type $SEQ_TYPE --strandedness $STRAND --rsem-genome-dir $RSEM_GENOME_DIR \
      --star-genome-dir $STAR_GENOME_DIR --zcat-flag --skip-txsort --skip-rsem
  else
    sbatch -n $CPU --mem ${MEMORY}G -o $LOG_DIR/RSEM_STAR_align_pipeline.STAR.${EXP_PREFIX}.log -J RSEM_STAR_align_pipeline_${EXP_PREFIX} \
      $PIPELINE_DIR/RSEM_STAR_align_pipeline.sh --thread $CPU --read1 $READ1 --read2 $READ2 --mem $MEMORY -o $ALIGN_DIR/$EXP_PREFIX \
      -p ${EXP_PREFIX} --seq-type $SEQ_TYPE --strandedness $STRAND --rsem-genome-dir $RSEM_GENOME_DIR \
      --star-genome-dir $STAR_GENOME_DIR --zcat-flag
  fi
}

function matrixToRun {
  ## $PARAM1 and $PARAM1 used for filter
  SAMPLE_MATRIX=$1
  PARAM1=$2
  PARAM2=$3
  test=$4
  i=1
  while IFS=$'\t' read -r -a matrixArr
  do
    test $i -eq 1 && ((i=i+1)) && continue
    ## GSE_accession
    param1="${matrixArr[3]}"
    ## GSM_accession
    param2="${matrixArr[4]}"
    exp_prefix="${matrixArr[9]}"
    library_layout="${matrixArr[14]}"
    strand="${matrixArr[15]}"
    read1="$FASTQ_DIR/${exp_prefix}_1.fastq.gz"
    read2="$FASTQ_DIR/${exp_prefix}_2.fastq.gz"
    if [[ $library_layout == "SE" ]]; then
      read2="none"
      if [[ ! -f $read1 ]]; then
        read1="$FASTQ_DIR/${exp_prefix}.fastq.gz"
      fi
    fi
    if [[ $test == "test" ]]; then
      RUNNING="runRSEM $exp_prefix $read1 $read2 $library_layout $strand test"
    else
      RUNNING="runRSEM $exp_prefix $read1 $read2 $library_layout $strand"
    fi
    if [[ $PARAM1 != "" ]]; then
      if [ $param1 == $PARAM1 ]; then
        if [[ $PARAM2 != "" ]]; then
          if [[ $param2 == $PARAM2 ]]; then
            $RUNNING
          fi
        else
          $RUNNING
        fi
      fi
    else
      $RUNNING
    fi
  done < $SAMPLE_MATRIX
}

CPU=8
MEMORY=80
########### alignment ##########
## normalB and B-ALL
#matrixToRun $SAMPLE_MATRIX GSE115656

## AML
#matrixToRun $SAMPLE_MATRIX GSE62852
#matrixToRun $SAMPLE_MATRIX GSE117994
#
### normal
#matrixToRun $SAMPLE_MATRIX GSE69239
#matrixToRun $SAMPLE_MATRIX GSE84022

#matrixToRun $SAMPLE_MATRIX GSE125112
#matrixToRun $SAMPLE_MATRIX GSE125437
matrixToRun $SAMPLE_MATRIX GSE54641
matrixToRun $SAMPLE_MATRIX GSE115988

## test
CPU=27
MEMORY=230
matrixToRun $SAMPLE_MATRIX GSE84022 GSM2225744
#matrixToRun $SAMPLE_MATRIX GSE62852 GSM1534507 test
#matrixToRun $SAMPLE_MATRIX GSE69239 GSM1695850 test
#matrixToRun $SAMPLE_MATRIX GSE84022 GSM2225743 test
#matrixToRun $SAMPLE_MATRIX GSE117994 GSM3316882 test
#matrixToRun $SAMPLE_MATRIX GSE125112 GSM3563631 test
#matrixToRun $SAMPLE_MATRIX GSE125437 GSM3574072 test
#matrixToRun $SAMPLE_MATRIX GSE54641 GSM1320481 test
#matrixToRun $SAMPLE_MATRIX GSE115988 GSM3204339 test
