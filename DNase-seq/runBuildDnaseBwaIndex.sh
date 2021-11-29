#!/bin/bash

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
# ref: https://github.com/ENCODE-DCC/dnase_pipeline/blob/master/dnanexus/dnase-index-bwa/resources/usr/bin/dnase_index_bwa.sh

function showHelp {
  echo -ne "usage: runChipAlignGenome.sh <options>\n"
  echo -e "options:
    -h | --help: show help information <bool>
    -b | --blacklist: the blacklist bed (optional) <str>
    -f | --fasta: genome fasta file <str>
    -g | --genome: the genome (hg38|hg19) <str>
    -m | --mappable: the mappable_only.starch file, download from ENCODE website (optional) <str>
    -o | --output: the output directory <str>
    -s | --script: the script directory ( bedops (bedmap,sort-bed,starch,starchcat,unstarch) ) <str>
    --skip: skip building the BWA index <bool>"
  exit 2;
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

TEMP=`getopt -o hb:f:g:m:o:, --long help,blacklist:,fasta:,genome:,mappable:,output:,script:, \
  --long skip \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

#echo "usage v1: dnase_index_bwa.sh <genome> <reference.fasta.gz> [<skip_indexing> <mappable_only.starch> [<blacklist.bed.gz>]]"

# Initialize variables:
BLACK_LIST="none"
GENOME=
GENOME_FASTA=
MAPPABLE_STARCH="none"
OUTPUT_DIR=
SCRIPT_DIR=
SKIP_INDEX=false

while true; do
  case "$1" in
    -h | --help ) showHelp; shift ;;
    -b | --blacklist ) BLACK_LIST="$2"; shift 2 ;;
    -f | --fasta ) GENOME_FASTA="$2"; shift 2 ;;
    -g | --genome ) GENOME="$2"; shift 2 ;;
    -m | --mappable ) MAPPABLE_STARCH="$2"; shift 2 ;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2 ;;
    -s | --script ) SCRIPT_DIR="$2"; shift 2 ;;
    --skip ) SKIP_INDEX=true; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done
## getopt end

## check required arguments

if [ -z "$OUTPUT_DIR" ]; then
  echo "--output: Wrong! output folder NOT found!"
  exit 2
fi

if [ -z "$GENOME_FASTA" ] || [ ! -f "$GENOME_FASTA" ]; then
  echo "--fasta: Wrong! input genome fasta NOT found!"
  exit 2
fi

if [ ! -f "$SCRIPT_DIR/bedops" ] && [ ! -f "$SCRIPT_DIR/unstarch" ]; then
  echo "--script: Wrong! script folder NOT found!"
  exit 2
fi

echo "Running pipleline with following parameters:"
echo "BLACK_LIST=$BLACK_LIST"
echo "GENOME=$GENOME"
echo "GENOME_FASTA=$GENOME_FASTA"
echo "MAPPABLE_STARCH=$MAPPABLE_STARCH"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo "SCRIPT_DIR=$SCRIPT_DIR"
echo "SKIP_INDEX=$SKIP_INDEX"

# add --script to the PATH
export PATH="$PATH:$SCRIPT_DIR/"

# mkdir output directory
if [ ! -d "$OUTPUT_DIR" ]; then
  echo "Create output direcotry!"
  mkdir -p $OUTPUT_DIR
fi

cd $OUTPUT_DIR

# copy mappable_only.starch and blacklist.bed to the output folder
cp $MAPPABLE_STARCH ./
cp $BLACK_LIST ./

genome=$GENOME               # Genome assembly (e.g. "GRCh38").
ref_fa_gz=$GENOME_FASTA            # Reference fasta file (e.g. GRCh38.fa.gz).  Will be uncompressed if not already.
skip_indexing=$SKIP_INDEX
mappable_only_starch=$(basename $MAPPABLE_STARCH)
blacklist_bed_gz=$(basename $BLACK_LIST)

index_file="${genome}_bwa_index.tgz"
echo "-- Index file will be: '$index_file'"

ref_fa=${ref_fa_gz%.gz}
if [ "$ref_fa" != "$ref_fa_gz" ]; then
    echo "-- Uncompressing reference"
    set -x
    cp $ref_fa_gz ./
    gunzip $ref_fa_gz
    mv $ref_fa ${genome}.fa
    set +x
fi

if [ "$ref_fa" != "$genome.fa" ]; then
    set -x
    cp $ref_fa ${genome}.fa
    set +x
fi

# Optionally create a mappable regions tar: faSize, sort-bed bedops starch unstarch extractCenterSites.sh
mappable_tar=""
if [ -f $mappable_only_starch ]; then
    mappable_tar="${genome}_hotspot2_v2.0_mappable.tgz"
    echo "-- Create hotspot2 mappable regions archive: ${mappable_tar}"
    echo "-- Generating chrom_sizes.bed from fasta file..."
    set -x
    faSize -detailed ${genome}.fa | awk '{printf "%s\t0\t%s\n",$1,$2}' | sort-bed - > chrom_sizes.bed
    set +x
    #cat $chrom_sizes | awk '{printf "%s\t0\t%s\n",$1,$2}' | sort-bed - > chrom_sizes.bed

    blacklist_bed=""
    if [ -f $blacklist_bed_gz ]; then
        blacklist_bed=${blacklist_bed_gz%.gz}
        if [ "$blacklist_bed" != "$blacklist_bed_gz" ]; then
            echo "-- Uncompressing blacklist..."
            set -x
            gunzip --stdout $blacklist_bed_gz | sort-bed - > $blacklist_bed
            set +x
        fi
        echo "-- Subtracting blacklist from mappable regions..."
        set -x
        bedops --difference $mappable_only_starch $blacklist_bed | starch - > mappable_target.starch
        set +x
    else
        echo "-- Using supplied mappable regions as mappable target file..."
        set -x
        cp $mappable_only_starch mappable_target.starch
        set +x
    fi

    echo "-- Creating centerSites file..."
    set -x
    extractCenterSites.sh -c chrom_sizes.bed -M mappable_target.starch -o center_sites.starch
    set +x

    echo "-- Archiving results"
    set -x
    cp $blacklist_bed ${genome}.blacklist.bed
    tar -czf $mappable_tar center_sites.starch mappable_target.starch chrom_sizes.bed $mappable_only_starch ${genome}.blacklist.bed
    set +x
fi

if [ "$skip_indexing" != "true" ]; then
    echo "-- Build index..."
    set -x
    bwa index -p $genome -a bwtsw ${genome}.fa
    set +x

    echo "-- tar and gzip index..."
    set -x
    tar -czf $index_file ${genome}.*
    set +x
        
    if [ "$ref_fa" != "${genome}.fa" ]; then
        set -x
        mv ${genome}.fa $ref_fa
        set +x
    fi
fi

echo "-- The results..."
if [ -f $index_file ]; then
    ls -l $index_file
fi
if [ -f $mappable_tar ]; then
    ls -l $mappable_tar
fi

#df -k .

