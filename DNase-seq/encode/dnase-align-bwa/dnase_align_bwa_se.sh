#!/bin/bash -e

if [ $# -ne 4 ]; then
    echo "usage v1: dnase_align_bwa_se.sh <index> <reads.fq.gz> <ncpus> <bam_root>"
    echo "Align single-end reads with bwa.  Is independent of DX and encodeD."
    echo "Requires bwa, edwBamStats, and samtools on path."
    exit -1; 
fi
index=$1   # BWA index archive including <ref_id>.fa (e.g. GRCh38.fa) and bwa index.
reads_fq_gz=$2 # fastq of of single-end reads, which will be trimmed resulting in "read1_trimmed.fq.gz"
ncpus=$3       # Number of cpus.
bam_root=$4    # root name for output bam (e.g. "out" will create "out.bam" and "out_flagstat.txt") 

echo "-- Expect to create '${bam_root}.bam'"

echo "-- Aligning with bwa..."
set -x
bwa aln -Y -l 32 -n 0.04 -k 2 -t $ncpus ${index} $reads_fq_gz > tmp.sai
bwa samse -n 10 ${index} tmp.sai $reads_fq_gz > tmp_se.sam
samtools view -Shb tmp_se.sam > tmp.bam
samtools sort --threads $ncpus -m 4G -o ${bam_root}.bam tmp.bam
samtools index ${bam_root}.bam
set +x

#rm -f tmp_se.sam tmp.bam tmp.sai

echo "-- Collect bam stats..."
set -x
samtools flagstat ${bam_root}.bam > ${bam_root}_flagstat.txt
edwBamStats ${bam_root}.bam ${bam_root}_edwBamStats.txt
set +x

echo "-- The results..."
ls -l ${bam_root}*
df -k .

