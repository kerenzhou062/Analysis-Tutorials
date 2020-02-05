#!/bin/bash
#SBATCH --job-name=DESeq2Run    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 1                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=10G                      # Amount of memory in GB
#SBATCH --time=120:10:00               # Time limit hrs:min:sec
#SBATCH --output=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/publicData/DESeq2/log/DESeq2Run.log               # Time limit hrs:min:sec

PIPELINE_DIR=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/home/public/pipeline/RNA-seq

BASE=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/publicData
LOG_DIR="$BASE/log"
ALIGNMENT_DIR="$BASE/alignment"
DESEQ2_DIR="$BASE/DESeq2"
SCRIPT="$BASE/scripts"
## genome
PUBLIC_BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome"
GENETYPE_FILE="$PUBLIC_BASE/annotation/hg38/geneType.hg38.txt"
#FULL_BED="$PUBLIC_BASE/annotation/hg38/v32/mix.gencodev32.tsnoRNA.anno.bed12"
FULL_BED="$PUBLIC_BASE/annotation/hg38/v32/gencode.v32.annotation.anno.bed12"
NCBI_GENE_INFO="$PUBLIC_BASE/annotation/hg38/hg38.ncbi.gene_info"

function DESeq2_Unbatch {
  COUNTS=$1
  SAMPLE=$2
  DESIGN=$3
  TREAT_ARR=$4
  CONTROL_ARR=$5
  PRERIX=$6
  OUTPUT_DIR=$7
  for ((i = 0; i < ${#TREAT_ARR[@]}; ++i));
  do
    TREAT=${TREAT_ARR[$i]}
    CONTROL=${CONTROL_ARR[$i]}
    # run RSEM_STAR_align_pipeline
    echo -e "\n$DESIGN: $TREAT vs $CONTROL, unbatch\n"
    DESeq2Gene.R --counts "$COUNTS" --sampleMtx "$SAMPLE" --design "$DESIGN" --control "$CONTROL" --treat "$TREAT" \
      --batchMethod 'none' --test 'Wald' --pval 1 --adjp 0.1 --shrink 'ashr' \
      --prefix "${PRERIX}_${TREAT}_${CONTROL}.unbatch" --output "$OUTPUT_DIR"
  done
}

function DESeq2_Batch {
  COUNTS=$1
  SAMPLE=$2
  DESIGN=$3
  TREAT_ARR=$4
  CONTROL_ARR=$5
  PRERIX=$6
  OUTPUT_DIR=$7
  echo "--ruvgCount $ruvgCount"
  echo "--ruvgLogfc $ruvgLogfc"
  echo "--ruvgPval $ruvgPval"
  for ((i = 0; i < ${#TREAT_ARR[@]}; ++i));
  do
    TREAT=${TREAT_ARR[$i]}
    CONTROL=${CONTROL_ARR[$i]}
    echo -e "\n$DESIGN: $TREAT vs $CONTROL, batch\n"
    # run RSEM_STAR_align_pipeline
    DESeq2Gene.R --counts "$COUNTS" --sampleMtx "$SAMPLE" --design "$DESIGN" --control "$CONTROL" --treat "$TREAT" \
      --batchMethod 'RUVg' --autoBatch --test 'Wald' --pval 1 --adjp 0.1 --shrink 'ashr' \
      --ruvgCount $ruvgCount --ruvgLogfc $ruvgLogfc --ruvgPval $ruvgPval \
      --prefix "${PRERIX}_${TREAT}_${CONTROL}.batch" --output "$OUTPUT_DIR"
  done
}

cd $DESEQ2_DIR

echo -e "running sample estimation\n"
sampleEvaluate.R --counts BCells.genes.matrix --sampleMtx BCells.all.sample.matrix --filter 1 --design1 'age' --design2 'cell' \
  --glmPca --normalize auto --poiHeatmap --prefix 'BCells_sampleEvaluate_age_cell' --output ./

sampleEvaluate.R --counts BCells.genes.matrix --sampleMtx BCells.all.sample.matrix --filter 1 --design1 'sex' --design2 'cell' \
  --glmPca --normalize auto --poiHeatmap --prefix 'BCells_sampleEvaluate_sex_cell' --output ./

## run DESeq2Time.R with all samples

$SCRIPT/DESeq2Time.R
