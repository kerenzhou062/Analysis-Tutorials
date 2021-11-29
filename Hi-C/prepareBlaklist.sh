#!/bin/bash
#SBATCH --job-name=prepareBlaklist    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 10                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=30G                      # Amount of memory in GB
#SBATCH --time=120:10:00               # Time limit hrs:min:sec
#SBATCH --output=prepareBlaklist.log               # Time limit hrs:min:sec

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
# ref:https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/DAC/STAR_RSEM_prep.sh

# https://bioinformaticsworkbook.org/dataAnalysis/GenomeAssembly/Hybrid/Scaffolding_with_HiC_Juicer.html

function showHelp {
  echo -ne "usage: sbatch prepareBlaklist.sh -t <thread_num> -o <log> --mem <30G> "
  echo -ne "prepareBlaklist.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -f | --fasta: genome fasta <str>
    -r | --rfasta: repeat sequence fasta <str>
    -t | --thread: # of cpus <int>
    -p | --prefix: output prefix <str>
    -o | --output: output directory <str>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o hf:r:t:, --long help,thread:,fasta:,rfasta:,gtf:, \
  --long rsem-genome-dir:,output:,spike-fa:,overhang:, \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREAD=10
FASTA=
REPEAT_FASTA=
PREFIX=
OUTPUT="./"

PREFIX="$(basename -- $FILE)"
PREFIX="${PREFIX%.fa}"
PREFIX="${PREFIX%.fasta}"

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -f | --fasta ) FASTA="$2"; shift 2 ;;
    -r | --rfasta ) REPEAT_FASTA="$2"; shift 2 ;;
    -t | --thread ) THREAD="$2"; shift 2 ;;
    -p | --prefix ) PREFIX="$2"; shift 2 ;;
    -o | --output ) OUTPUT="$2"; shift 2 ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

## getopt end

## check required arguments

if [ -z $FASTA ]; then
  echo "-f|--fasta: Wrong! Genome fasta NOT found!"
  exit 2
fi

if [ -z $REPEAT_FASTA ]; then
  echo "-rf|--rfasta: Wrong! Repeat sequence fasta NOT found!"
  exit 2
fi

if [ -d $OUTPUT ]; then
  echo "-o | --output: Warning! Output folder NOT found!"
  mkdir -p $OUTPUT
fi

echo "Running pipleline with following parameters:"
echo "THREAD=$THREAD"
echo "FASTA=$FASTA"
echo "PREFIX=$PREFIX"
echo "OUTPUT=$OUTPUT"
echo ""

OUTPUT_FILE_PREFIX="$OUTPUT/$PREFIX"
DB_FILE="$OUTPUT_FILE_PREFIX.DB"
REPEATS_FASTA="$OUTPUT_FILE_PREFIX.repeats.fasta"
BLAST_OUT="$OUTPUT_FILE_PREFIX.repeats2Genome.blastout"
REPEAT_BED="$OUTPUT_FILE_PREFIX.repeats2Genome.blastout"

echo "Making blastdb..."
echo "makeblastdb -in $FASTA -dbtype nucl -out $DB_FILE"

makeblastdb -in $FASTA -dbtype nucl -out $DB_FILE

echo "Blasting blastdb..."
echo "makeblastdb -in $FASTA -dbtype nucl -out $DB_FILE"
blastn -db $DB_FILE -dust no -num_threads 16 -outfmt 6 -query $REPEATS_FASTA -evalue 100 -num_alignments 100000 -out $BLAST_OUT

cat $BLAST_OUT |awk '{if($10>$9){print $2,$9,$10} else {print $2,$10,$9}}' |tr " " "\t" |sort -k1,1V -k2,3n | \
  uniq | bedtools merge -d 1000 -i stdin > $REPEAT_BED
