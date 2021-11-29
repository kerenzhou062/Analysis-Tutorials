#!/bin/sh
BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/brandontan/ythdc1/rna_seq/inhouse_data"

EXP_DIR="$BASE/DESeq2/nuc_cyto"
# GENOME
PUBLIC_BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome"
GENETYPE_FILE="$PUBLIC_BASE/annotation/hg38/geneType.hg38.txt"
NCBI_GENE_INFO="$PUBLIC_BASE/annotation/hg38/hg38.ncbi.gene_info"
FULL_BED="$PUBLIC_BASE/annotation/hg38/v33/mix.gencodev33.tsnoRNA.anno.bed12"

if [[ ! -d "$EXP_DIR" ]]; then
  mkdir -p $EXP_DIR
fi
rm -rf $EXP_DIR/batch
rm -rf $EXP_DIR/unbatch
#mkdir $EXP_DIR
cd $EXP_DIR

##get counts of genes
buildExpMatrix.py -input $BASE/alignment -identity genes -grepKept "nuc|cyto" -abundance counts \
  -output rna.gene.counts.matrix

##get TPM of isoforms
buildExpMatrix.py -input $BASE/alignment -identity genes -grepKept "nuc|cyto" -abundance TPM \
  -output tmp.TPM.matrix

appendGeneAnno.py -input tmp.TPM.matrix -idType gene -idCol 0 -inCol 0 -anno $FULL_BED \
  -geneClassFile $GENETYPE_FILE -ncbiGeneInfo $NCBI_GENE_INFO -output rna.genes.TPM.matrix

##get expression of isoforms
buildExpMatrix.py -input $BASE/alignment -identity isoforms -grepKept "nuc|cyto" -abundance TPM \
  -output tmp.TPM.matrix

appendGeneAnno.py -input tmp.TPM.matrix -idType tx -idCol 0 -inCol 0 -anno $FULL_BED \
  -geneClassFile $GENETYPE_FILE -ncbiGeneInfo $NCBI_GENE_INFO -output rna.isoforms.TPM.matrix

rm -f tmp.TPM.matrix
# evaluate sample distances

sampleEvaluate.R --counts rna.gene.counts.matrix --sampleMtx rna.sample.nuc_cyto.matrix \
  --filter 1 --design1 'condition' --design2 'replicate' \
  --normalize auto --poiHeatmap --prefix 'ythdc1.nuc_cyto' --output ./


function AppendGeneName {
  prefix=$1
  directory=$2
  output="$prefix.unbatch.DESeq2.txt"
  sigoutput="$prefix.unbatch.sig.DESeq2.txt"

  rename $directory/$prefix $directory/tmp $directory/${prefix}*.txt
  appendGeneAnno.py -input $directory/tmp.DESeq2.txt -idCol 0 -inCol 0 -anno $FULL_BED \
    -geneClassFile $GENETYPE_FILE -ncbiGeneInfo $NCBI_GENE_INFO -output $directory/$output
  
  appendGeneAnno.py -input $directory/tmp.sig.DESeq2.txt -idCol 0 -inCol 0 -anno $FULL_BED \
    -geneClassFile $GENETYPE_FILE -ncbiGeneInfo $NCBI_GENE_INFO -output $directory/$sigoutput
  ## delete temp files
  rm -f $directory/tmp*
}
## differential gene expression analysis
function RunDESeq2 {
  matrix=$1
  counts=$2
  control=$3
  treat=$4
  design=$5
  prefix=$6
  directory=$7
  batch=$8
  args=""
  output="$prefix.unbatch.DESeq2.txt"
  sigoutput="$prefix.unbatch.sig.DESeq2.txt"
  ## if run with batch effect correction
  if [[ "$batch" == "true" ]]; then
    args=" --batchMethod spikeins "
    output="$prefix.batch.DESeq2.txt"
    sigoutput="$prefix.batch.sig.DESeq2.txt"
  fi

  DESeq2Gene.R --sampleMtx $matrix --counts $counts \
    --control $control --treat $treat --design $design -f 5 $args \
    --test 'Wald' --pval 1 --adjp 0.1 --shrink 'ashr' --prefix $prefix --output $directory
  AppendGeneName $prefix $output_dir
}

echo "perform without batch corrections"
output_dir="unbatch"
mkdir -p $output_dir

RunDESeq2 rna.sample.nuc_cyto.matrix rna.gene.counts.matrix "sgNS_cyto" "sgA_cyto" "condition" "ythdc1_cyto_sgA_sgNS" $output_dir false
RunDESeq2 rna.sample.nuc_cyto.matrix rna.gene.counts.matrix "sgNS_cyto" "sgB_cyto" "condition" "ythdc1_cyto_sgB_sgNS" $output_dir false

RunDESeq2 rna.sample.nuc_cyto.matrix rna.gene.counts.matrix "sgNS_nuc" "sgA_cyto" "condition" "ythdc1_nuc_sgA_sgNS" $output_dir false
RunDESeq2 rna.sample.nuc_cyto.matrix rna.gene.counts.matrix "sgNS_nuc" "sgB_cyto" "condition" "ythdc1_nuc_sgB_sgNS" $output_dir false

prefix="sgA_nucVScyto_vs_sgNS_nucVScyto"
DESeq2Gene.R --sampleMtx rna.sample.nuc_cyto.matrix --counts rna.gene.counts.matrix \
  --formula "phenotype + fraction + phenotype:fraction" --relevel "phenotype:sgNS,fraction:cyto" \
  --name "phenotypesgA.fractionnuc" -f 5 \
  --test 'LRT' --pval 1 --adjp 0.1 --shrink 'none' --prefix $prefix --output $output_dir
AppendGeneName $prefix $output_dir

prefix="sgB_nucVScyto_vs_sgNS_nucVScyto"
DESeq2Gene.R --sampleMtx rna.sample.nuc_cyto.matrix --counts rna.gene.counts.matrix \
  --formula "phenotype + fraction + phenotype:fraction" --relevel "phenotype:sgNS,fraction:cyto" \
  --name "phenotypesgB.fractionnuc" -f 5 \
  --test 'LRT' --pval 1 --adjp 0.1 --shrink 'none' --prefix "sgB_nucVScyto_vs_sgNS_nucVScyto" --output $output_dir
AppendGeneName $prefix $output_dir

echo "perform with batch corrections"
output_dir="batch"
mkdir -p $output_dir

RunDESeq2 rna.sample.nuc_cyto.matrix rna.gene.counts.matrix "sgNS_cyto" "sgA_cyto" "condition" "ythdc1_cyto_sgA_sgNS" $output_dir true
RunDESeq2 rna.sample.nuc_cyto.matrix rna.gene.counts.matrix "sgNS_cyto" "sgB_cyto" "condition" "ythdc1_cyto_sgB_sgNS" $output_dir true

RunDESeq2 rna.sample.nuc_cyto.matrix rna.gene.counts.matrix "sgNS_nuc" "sgA_cyto" "condition" "ythdc1_nuc_sgA_sgNS" $output_dir true
RunDESeq2 rna.sample.nuc_cyto.matrix rna.gene.counts.matrix "sgNS_nuc" "sgB_cyto" "condition" "ythdc1_nuc_sgB_sgNS" $output_dir true

prefix="sgA_nucVScyto_vs_sgNS_nucVScyto"
DESeq2Gene.R --sampleMtx rna.sample.nuc_cyto.matrix --counts rna.gene.counts.matrix --batchMethod "spikeins" \
  --formula "phenotype + fraction + phenotype:fraction" --relevel "phenotype:sgNS,fraction:cyto" \
  --name "phenotypesgA.fractionnuc" -f 5  \
  --test 'LRT' --pval 1 --adjp 0.1 --shrink 'none' --prefix "sgA_nucVScyto_vs_sgNS_nucVScyto" --output $output_dir
AppendGeneName $prefix $output_dir

prefix="sgB_nucVScyto_vs_sgNS_nucVScyto"
DESeq2Gene.R --sampleMtx rna.sample.nuc_cyto.matrix --counts rna.gene.counts.matrix --batchMethod "spikeins" \
  --formula "phenotype + fraction + phenotype:fraction" --relevel "phenotype:sgNS,fraction:cyto" \
  --name "phenotypesgB.fractionnuc" -f 5  \
  --test 'LRT' --pval 1 --adjp 0.1 --shrink 'none' --prefix "sgB_nucVScyto_vs_sgNS_nucVScyto" --output $output_dir
AppendGeneName $prefix $output_dir
