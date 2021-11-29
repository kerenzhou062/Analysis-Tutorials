#!/bin/bash
BASE=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/publicData
LOG_DIR="$BASE/log"
ALIGNMENT_DIR="$BASE/alignment"
EXP_DIR="$BASE/expression"
SCRIPT="$BASE/scripts"
## genome
PUBLIC_BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome"
GENETYPE_FILE="$PUBLIC_BASE/annotation/hg38/geneType.hg38.txt"
FULL_BED="$PUBLIC_BASE/annotation/hg38/v32/mix.gencodev32.tsnoRNA.anno.bed12"

cd $EXP_DIR


buildExpMatrix.py -input $ALIGNMENT_DIR -output exp.tmp -abundance TPM
appendGeneAnno.py -input exp.tmp -output exression.genes.tpm.matrix \
  -idCol 0 -inCol 0 -anno $FULL_BED -onlyName

geneExpToLong.py -input BCells.genes.tpm.matrix -matrix BCells.all.sample.matrix -output BCells.genes.tpm.long.matrix

drawExpBoxplot.R --input GSE28497_rma.long.matrix.gz --gzip --output ./ \
  --targetCol gene_name --target "TET1,TET2,TET3,FTO,METTL3,WTAP,YTHDF2,YTHDC1,IGF2BP1,IGF2BP2,IGF2BP3" --prefix RMA \
  --xaxis cancerType --yaxis expression --refgroup NBM 

