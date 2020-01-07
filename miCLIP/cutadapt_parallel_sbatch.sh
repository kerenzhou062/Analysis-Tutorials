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
THREAD=$n
FASTQ_DIR=`realpath ./fastq`
CUTADAPT_DIR=`realpath ./cutadapt`
# parallel read preprocessing
if [[ ! -d $CUTADAPT_DIR ]]; then
  mkdir $CUTADAPT_DIR
fi

rm -rf tmpfifo
mkfifo tmpfifo
exec 9<>tmpfifo

for((n=1;n<=${THREAD};n++))
do
  echo >&9
done

startTime=$(date "+%s")

for i in $FASTQ_DIR/*.fastq;
do
  read -u9
  {
    PREFIX=${i%%.fastq}
    PREFIX=${PREFIX##*/}
    cutadapt -f fastq --times 2 -e 0 -O 1 \
      --quality-cutoff 6 -m 13 -a AGATCGGAAGAGCACACGTCT -a AAAAAAAAAA \
      -o ${CUTADAPT_DIR}/${PREFIX}.trim.fastq $i \
      > ${CUTADAPT_DIR}/${PREFIX}.cutadpt.log 2>&1
    echo >&9  
  } &
done

wait

endTime=$(date "+%s")
echo "Run Time:" $(($endTime - $startTime))"s"
echo "cutadpt done!"

exec 9>&-
exec 9<&-
rm -rf tmpfifo
## parallel end
