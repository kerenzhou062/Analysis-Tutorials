#!/bin/bash
#SBATCH --job-name=cutadapt    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=50G                     # Amount of memory in GB
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/hlhuang/SETD2/miCLIP/cutadapt.log
#serial_test_%j.logserial_test_%j.log

# basic variables
BASE=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/hlhuang/SETD2/miCLIP
THREAD=$n
FASTQ_DIR="$BASE/fastq"
CUTADAPT_DIR="$BASE/cutadapt"
# parallel read preprocessing
if [[ ! -d $CUTADAPT_DIR ]]; then
    mkdir $CUTADAPT_DIR
fi

#for i in $FASTQ_DIR/*.fastq;
#do
#    PREFIX=${i%%.fastq}
#    PREFIX=${PREFIX##*/}
#    srun -o ${CUTADAPT_DIR}/${PREFIX}.cutadpt.log --mem 10G \
#      cutadapt -f fastq --times 2 -e 0 -O 1 \
#        --quality-cutoff 6 -m 18 -a AGATCGGAAGAGCACACGTCT -a AAAAAAAAAA \
#        -o ${CUTADAPT_DIR}/${PREFIX}.trim.fastq $i &
#done

for i in $FASTQ_DIR/*.fastq;
do
    PREFIX=${i%%.fastq}
    PREFIX=${PREFIX##*/}
    srun -o ${CUTADAPT_DIR}/${PREFIX}.cutadpt.log --mem 10G \
      cutadapt -f fastq --times 2 -e 0 -O 3 \
        --quality-cutoff 6 -m 13 -a AGATCGGAAGAGCACACGTCT -a AAAAAAAAAA \
        -o ${CUTADAPT_DIR}/${PREFIX}.trim.fastq $i &
done
