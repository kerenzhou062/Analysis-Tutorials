#!/bin/sh

cd /net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/publicData/fastq
for i in *.fastq;
do
  srun --mem 5G gzip $i &
done
