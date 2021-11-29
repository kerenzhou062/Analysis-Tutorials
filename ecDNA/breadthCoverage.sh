#!/bin/bash
#SBATCH --job-name=breadthCoverage    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 1                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=10G                      # Amount of memory in GB
#SBATCH --time=60:00:00               # Time limit hrs:min:sec
#SBATCH --output=breadthCoverage.log               # Time limit hrs:min:sec

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.

function showHelp {
  echo -ne "usage: sbatch breadthCoverage.sh -n 1 -o <log> --mem <60G> "
  echo -ne "breadthCoverage.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -b | --bam: sorted bam <str>
    -c | --coverage: minimum coverage depth [5] <float>
    -g | --gsize: genome size file <str>
    -o | --output: output file <str>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o hb:c:g:o, --long coverage:,bam:,gsize:,output:, \
  -- "$@"`


if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Initialize variables:
COVERAGE=5
BAM_FILE=60
GSIZE_FILE=
OUTPUT="result.txt"

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -b | --bam ) BAM_FILE="$2"; shift 2 ;;
    -c | --coverage ) COVERAGE="$2"; shift 2 ;;
    -g | --gsize ) GSIZE_FILE="$2"; shift 2 ;;
    -o | --output ) OUTPUT="$2"; shift 2 ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

## getopt end

## check required arguments
if [ -z $GSIZE_FILE ] && [ ! -f $GSIZE_FILE ] ; then
  echo "-g|--gsize: Wrong! genome size file NOT found!"
  exit 2
fi

if [ -z $BAM_FILE ] && [ ! -f $BAM_FILE ] ; then
  echo "-b|--bam: Wrong! sorted bam file NOT found!"
  exit 2
fi

if [ ! -f "$BAM_FILE.bai" ] ; then
  echo "-b|--bam: Wrong! sorted bam index NOT found!"
  exit 2
fi

echo "Running pipleline with following parameters:"
echo "BAM_FILE=$BAM_FILE"
echo "COVERAGE=$COVERAGE"
echo "GSIZE_FILE=$GSIZE_FILE"
echo "OUTPUT=$OUTPUT"
echo ""

base_num=`samtools mpileup --count-orphans $BAM_FILE | awk -v X="${COVERAGE}" '{if($4>=X) print}' | wc -l`
ref_length=`awk '{ FS = "\t"} ; BEGIN{L=0}; {L=L+$2}; END{print L}' $GSIZE_FILE`

echo "$base_num $ref_length" | awk '{p=$1/$2;print "breadthCoverage\t"p*100}' > $OUTPUT
