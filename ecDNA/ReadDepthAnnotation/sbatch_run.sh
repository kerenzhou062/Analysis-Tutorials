#!/bin/bash
#arguments

readLength=100
fastaDir="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome/fasta/hg38/chroms"
outdir="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome/programData/readDepth/annotations/100bp"

cd /net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome/programData/readDepth

chrom_arr=( "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10")
chrom_arr+=( "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20")
chrom_arr+=( "chr21" "chr22" "chrM" "chrX" "chrY" )

for ((i = 0; i < ${#chrom_arr[@]}; ++i));
do
  chr=${chrom_arr[$i]}
  sbatch -J ${chr}_readDepth -o ${chr}_readDepth.log ./runEachChr.sh $chr $readLength $fastaDir $outdir
done

readLength=150
fastaDir="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome/fasta/hg38/chroms"
outdir="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome/programData/readDepth/annotations/150bp"

cd /net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome/programData/readDepth

chrom_arr=( "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10")
chrom_arr+=( "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20")
chrom_arr+=( "chr21" "chr22" "chrM" "chrX" "chrY" )

for ((i = 0; i < ${#chrom_arr[@]}; ++i));
do
  chr=${chrom_arr[$i]}
  sbatch -J ${chr}_readDepth -o ${chr}_readDepth.log ./runEachChr.sh $chr $readLength $fastaDir $outdir
done
