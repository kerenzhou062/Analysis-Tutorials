#!/bin/sh
BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/publicData"
FASTQ_DIR="$BASE/fastq"

cd $FASTQ_DIR

function matrixToRun {
  ## $PARAM1 and $PARAM1 used for filter
  SAMPLE_MATRIX=$1
  i=1
  while IFS=$'\t' read -r -a matrixArr
  do
    test $i -eq 1 && ((i=i+1)) && continue
    srr="${matrixArr[0]}"
    exp_prefix="${matrixArr[9]}"
    rename "$srr" "${exp_prefix}" ${srr}*.fastq
  done < $SAMPLE_MATRIX
}

matrixToRun $BASE/sample.multipleRun.matrix
