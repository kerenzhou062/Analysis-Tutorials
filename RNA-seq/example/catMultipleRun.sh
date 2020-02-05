#!/bin/sh
BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/publicData"
FASTQ_DIR="$BASE/fastq"

cd $FASTQ_DIR

srun --mem 4G cat RS411_24h_rep1_run*.fastq > RS411_24h_rep1.fastq &
srun --mem 4G cat RS411_24h_rep2_run*.fastq > RS411_24h_rep2.fastq &
srun --mem 4G cat RS411_24h_rep3_run*.fastq > RS411_24h_rep3.fastq &
srun --mem 4G cat RS411_24h_rep4_run*.fastq > RS411_24h_rep4.fastq &
srun --mem 4G cat RS411_24h_rep5_run*.fastq > RS411_24h_rep5.fastq &
srun --mem 4G cat RS411_24h_rep6_run*.fastq > RS411_24h_rep6.fastq &
srun --mem 4G cat RS411_96h_rep1_run*.fastq > RS411_96h_rep1.fastq &
srun --mem 4G cat RS411_96h_rep2_run*.fastq > RS411_96h_rep2.fastq &
srun --mem 4G cat RS411_96h_rep3_run*.fastq > RS411_96h_rep3.fastq &
srun --mem 4G cat SEM_24h_rep1_run*.fastq > SEM_24h_rep1.fastq &
srun --mem 4G cat SEM_24h_rep2_run*.fastq > SEM_24h_rep2.fastq &
srun --mem 4G cat SEM_24h_rep3_run*.fastq > SEM_24h_rep3.fastq &
srun --mem 4G cat SEM_24h_rep4_run*.fastq > SEM_24h_rep4.fastq &
srun --mem 4G cat SEM_24h_rep5_run*.fastq > SEM_24h_rep5.fastq &
srun --mem 4G cat SEM_24h_rep6_run*.fastq > SEM_24h_rep6.fastq &
srun --mem 4G cat SEM_96h_rep1_run*.fastq > SEM_96h_rep1.fastq &
srun --mem 4G cat SEM_96h_rep2_run*.fastq > SEM_96h_rep2.fastq &
srun --mem 4G cat SEM_96h_rep3_run*.fastq > SEM_96h_rep3.fastq &
srun --mem 4G cat KOPN8_24h_rep1_run*.fastq > KOPN8_24h_rep1.fastq &
srun --mem 4G cat KOPN8_24h_rep2_run*.fastq > KOPN8_24h_rep2.fastq &
srun --mem 4G cat KOPN8_24h_rep3_run*.fastq > KOPN8_24h_rep3.fastq &
srun --mem 4G cat KOPN8_24h_rep4_run*.fastq > KOPN8_24h_rep4.fastq &
srun --mem 4G cat KOPN8_24h_rep5_run*.fastq > KOPN8_24h_rep5.fastq &
srun --mem 4G cat KOPN8_24h_rep6_run*.fastq > KOPN8_24h_rep6.fastq &
srun --mem 4G cat KOPN8_96h_rep1_run*.fastq > KOPN8_96h_rep1.fastq &
srun --mem 4G cat KOPN8_96h_rep2_run*.fastq > KOPN8_96h_rep2.fastq &
srun --mem 4G cat KOPN8_96h_rep3_run*.fastq > KOPN8_96h_rep3.fastq &
