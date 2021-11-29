#!/bin/bash
#SBATCH --job-name=createCustomAnnotations    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 6                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=60G                      # Amount of memory in GB
#SBATCH --time=240:10:00               # Time limit hrs:min:sec
#SBATCH --output=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/krzhou/ecDNA/wgs_data/publicData/ReadDepth/createCustomAnnotations.log               # Time limit hrs:min:sec

#arguments
chr=$1
readLength=$2
fastaDir=$3
outdir=$4
windowSize=100
BWA_INDEX="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome/index/bwa/hg38.fa"
GENOME_SIZE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome/size/hg38.chrom.sizes"
SCRIPT="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome/programData/readDepth"

#generate all possible reads of the given size from that chromsome
perl $SCRIPT/allSeq.pl $fastaDir/$chr.fa $readLength > $outdir/$chr.fastq
# align those reads using bwa
bwa aln -t 6 $BWA_INDEX $outdir/$chr.fastq >$outdir/$chr.aln.sai
bwa samse -n 1 $BWA_INDEX $outdir/$chr.aln.sai $outdir/$chr.fastq >$outdir/$chr.sam
#clean up
rm -f $outdir/$chr.fastq $outdir/$chr.aln.sai

#convert the sam file to a mapability and gc-content tracks 
perl $SCRIPT/mapAndGc.pl $outdir/$chr.sam $GENOME_SIZE $windowSize $chr $readLength $outdir

mv $outdir/$chr.map $outdir/$chr.dat

gzip $outdir/$chr.gc
gzip $outdir/$chr.dat

mv $outdir/$chr.gc.gz $outdir/gcWinds
mv $outdir/$chr.dat.gz $outdir/mapability

rm -f $outdir/$chr.sam
