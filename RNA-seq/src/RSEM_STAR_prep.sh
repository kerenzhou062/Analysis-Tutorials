#!/bin/bash
#SBATCH --job-name=RSEM_STAR_prep    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 10                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=30G                      # Amount of memory in GB
#SBATCH --time=72:10:00               # Time limit hrs:min:sec
#SBATCH --output=RSEM_STAR_prep.log               # Time limit hrs:min:sec

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
# ref:https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/DAC/STAR_RSEM_prep.sh

# output: all in the working directory, fixed names
# Aligned.sortedByCoord.out.bam                 # alignments, standard sorted BAM, agreed upon formatting
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Quant.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
# Quant.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
# Quant.pdf                                     # RSEM diagnostic plots
# Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data
# Signal.{Unique,UniqueMultiple}.unstranded.bw  # 2 bigWig files for unstranded data

function showHelp {
  echo -ne "usage: sbatch RSEM_STAR_prep.sh -n <thread_num> -o <log> --mem <30G> "
  echo -ne "RSEM_STAR_prep.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -f | --fasta: genome fasta <str>
    -t | --thread: # of cpus <int>
    --overhang: --sjdbOverhang (100) <float>
    --gtf: reference gtf (contained spike-ins) <int>
    --spike-fa: fasta of spike-ins <int>
    --rsem-genome-dir: RSEM genome directory <str>
    --star-genome-dir: STAR genome directory <str>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o hf:t: --long help,thread:,fasta:,gtf: \
  --long rsem-genome-dir:,star-genome-dir:,spike-fa:,overhang: \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
THREAD=10
GENOME_FASTA=
STAR_GENOME_DIR=
RSEM_GENOME_DIR=
GTF=
OVERHANG=100
SPIKE_FASTA=

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -f | --fasta ) GENOME_FASTA="$2"; shift 2 ;;
    -t | --thread ) THREAD="$2"; shift 2 ;;
    --rsem-genome-dir ) RSEM_GENOME_DIR="$2"; shift 2 ;;
    --star-genome-dir ) STAR_GENOME_DIR="$2"; shift 2 ;;
    --overhang ) OVERHANG="$2"; shift 2 ;;
    --gtf ) GTF="$2"; shift 2 ;;
    --spike-fa ) SPIKE_FASTA="$2"; shift 2 ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

## getopt end

## check required arguments
if [ -z $SPIKE_FASTA ]; then
  echo "--SPIKE_FASTA: Wrong! Spike-in FASTA NOT found!"
  exit 2
fi

if [ -z $GENOME_FASTA ]; then
  echo "-f|--fasta: Wrong! Genome fasta NOT found!"
  exit 2
fi

if [ -z $RSEM_GENOME_DIR ]; then
  echo "--rsem-genome-dir: Wrong! RSEM genome folder NOT found!"
  exit 2
fi

if [ -z $STAR_GENOME_DIR ]; then
  echo "--star-genome-dir: Wrong! STAR genome folder NOT found!"
  exit 2
fi

if [ -z ${GTF} ]; then
  echo "--gtf: Wrong! gtf file NOT found!"
  exit 2
fi

echo "Running pipleline with following parameters:"
echo "THREAD=$THREAD"
echo "GENOME_FASTA=$GENOME_FASTA"
echo "STAR_GENOME_DIR=$STAR_GENOME_DIR"
echo "RSEM_GENOME_DIR=$RSEM_GENOME_DIR"
echo "GTF=${GTF}"
echo "OVERHANG=$OVERHANG"
echo "SPIKE_FASTA=$SPIKE_FASTA"
echo ""

# example
# ./STAR_RSEM_prep.sh  /path/to/STARgenome  /path/to/RSEMgenome  male.hg19.fa  spikes.fixed.fasta gencode.v19.annotation_tRNA_spikeins.gtf

# RSEM genome
if [[ ! -d $RSEM_GENOME_DIR ]]; then
  mkdir -p $RSEM_GENOME_DIR
else
  rm -rf $RSEM_GENOME_DIR/*
fi

# STAR genome
if [[ ! -d $STAR_GENOME_DIR ]]; then
  mkdir -p $STAR_GENOME_DIR
else
  rm -rf $STAR_GENOME_DIR/*
fi

### the command below is for RSEM >=1.2.19
### note, that for RSEM < 1.2.19, --no-polyA should be added

RSEMcommand="rsem-prepare-reference --gtf ${GTF} ${GENOME_FASTA}","${SPIKE_FASTA} $RSEM_GENOME_DIR/RSEMref"
echo $RSEMcommand
$RSEMcommand


# STAR genome
mkdir $STAR_GENOME_DIR
STARcommand="STAR --runThreadN ${THREAD} --runMode genomeGenerate --genomeDir ${STAR_GENOME_DIR} --genomeFastaFiles ${GENOME_FASTA} \
  ${SPIKE_FASTA} --sjdbGTFfile ${GTF} --sjdbOverhang ${OVERHANG} --outFileNamePrefix $STAR_GENOME_DIR"
echo $STARcommand
$STARcommand
