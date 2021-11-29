#!/bin/sh
BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/publicData"
FASTQ_DIR="$BASE/fastq"

cd $FASTQ_DIR

for i in *.sra;
do
  srun --mem 5G fastq-dump --split-files $i &
done
