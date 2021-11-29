#!/bin/sh
#SBATCH --job-name=geneEnrich_ChIPRep2_clusterProfiler4    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 5                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=200G                      # Amount of memory in GB
#SBATCH --time=120:10:00               # Time limit hrs:min:sec
#SBATCH --output=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/inhouse_data/log/geneEnrich_ChIPRep2_clusterProfiler4.PDX2.sig.log               # Time limit hrs:min:sec

function MultiProcesser {
  commandArr="${1}"
  # parallel start
  THREAD=5
  tmpfifo="/tmp/$RANDOM.fifo"
  mkfifo tmpfifo
  exec 9<>tmpfifo
  rm -rf tmpfifo

  for((n=1;n<=${THREAD};n++))
  do
    echo >&9
  done

  for command in "${commandArr[@]}"
  do
    read -u9
    {
      echo $command
      eval "$command"
      echo >&9
    } &
  done
  wait
  exec 9>&-
  exec 9<&-
  # parallel end
}

function runGSEA {
  input=$1
  result_dir=$2
  item=$3
  prefix=$4
  maxGSSize=$5
  method=$6
  mkdir -p $result_dir
  declare -a categoryArr=( "H" "C2" "C3" "C4" "C5" "C6" "C7" )
  declare -a commandArr=()
  for category in "${categoryArr[@]}"
  do
    command="runClusterProfiler4.R --category $category --input $input --minGSSize 5 --maxGSSize $maxGSSize \
     --organism human --pval 0.05 --type GSEA --output $result_dir/$category \
     --item $item --prefix "$prefix" --by $method --nperm 1000 > $result_dir/$prefix.$category.log 2>&1"
    commandArr+=("$command")
  done
  MultiProcesser "$commandArr"
}

function runGO {
  input=$1
  result_dir=$2
  item=$3
  prefix=$4
  maxGSSize=$5
  mkdir -p $result_dir
  declare -a gotypeArr=( "BP" "MF" "CC" "ALL")
  declare -a commandArr=()
  for gotype in "${gotypeArr[@]}"
  do
    command="runClusterProfiler4.R --input $input \
       --organism human --pval 0.05 --qval 0.1 --type GO --output $result_dir --maxGSSize $maxGSSize \
       --gotype $gotype --max 30 --noCnet --noHeatmap --item $item --prefix "$prefix.$gotype" \
       > $result_dir/$prefix.$gotype.log 2>&1"
    commandArr+=("$command")
  done
  MultiProcesser "$commandArr"
}

function runOther {
  input=$1
  root_dir=$2
  item=$3
  prefix=$4
  maxGSSize=$5
  declare -a typeArr=( "KEGG" "wiki" "reactome" "disease")
  declare -a commandArr=()
  for type in "${typeArr[@]}"
  do
    result_dir="${root_dir}/$type/${maxGSSize}_size"
    if [[ ! -d "$result_dir" ]]; then
      mkdir -p $result_dir
    fi
    command="runClusterProfiler4.R --input $input \
      --organism human --pval 0.05 --qval 0.1 --type $type --output $result_dir --maxGSSize $maxGSSize \
      --max 30 --item $item --prefix "$prefix.$type" > $result_dir/$prefix.$type.log 2>&1"
    commandArr+=("$command")
  done
  MultiProcesser "$commandArr"
}

DESeq2_DIR="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/inhouse_data/DESeq2/PDX2/unbatch"
ENRICH_DIR="$DESeq2_DIR/geneEnrich_ChIPRep2_sig_clusterProfiler4"
mkdir -p $ENRICH_DIR
cd $ENRICH_DIR

## for GO and KEGG

## down regulate
awk 'BEGIN{FS="\t";OFS="\t"}
ARGIND==1{
  geneArr[$1]=$1
}
ARGIND==2{
  if (FNR >1) {
    if ($8 < 0 && $1 in geneArr ) {
      print $1, $8
    }
  }
}' $DESeq2_DIR/../PDX2_TET1_rep2.ChipInReps.gene.txt $DESeq2_DIR/PDX2_TET1_g1g2g3_gcontrol.unbatch.sig.DESeq2.txt > PDX2_TET1_g1g2g3.sig.down.txt

## up regulated
awk 'BEGIN{FS="\t";OFS="\t"}
ARGIND==1{
  geneArr[$1]=$1
}
ARGIND==2{
  if (FNR >1) {
    if ($8 > 0 && $1 in geneArr ) {
      print $1, $8
    }
  }
}' $DESeq2_DIR/../PDX2_TET1_rep2.ChipInReps.gene.txt $DESeq2_DIR/PDX2_TET1_g1g2g3_gcontrol.unbatch.sig.DESeq2.txt > PDX2_TET1_g1g2g3.sig.up.txt

## significant
awk 'BEGIN{FS="\t";OFS="\t"}
ARGIND==1{
  geneArr[$1]=$1
}
ARGIND==2{
  if (FNR >1) {
    if ($1 in geneArr ) {
      print $1, $8
    }
  }
}' $DESeq2_DIR/../PDX2_TET1_rep2.ChipInReps.gene.txt $DESeq2_DIR/PDX2_TET1_g1g2g3_gcontrol.unbatch.sig.DESeq2.txt > PDX2_TET1_g1g2g3.sig.txt


RESULT_DIR="GSEA/PDX2_TET1_g1g2g3_sig/fgsea"
runGSEA PDX2_TET1_g1g2g3.sig.txt $RESULT_DIR ENSEMBL "PDX2_TET1_g1g2g3_sig" 500 fgsea
RESULT_DIR="GSEA/PDX2_TET1_g1g2g3_sig/DOSE"
runGSEA PDX2_TET1_g1g2g3.sig.txt $RESULT_DIR ENSEMBL "PDX2_TET1_g1g2g3_sig" 500 DOSE


RESULT_DIR="GO/500_size"
runGO PDX2_TET1_g1g2g3.sig.down.txt $RESULT_DIR ENSEMBL "PDX2_TET1_g1g2g3_sig_down" 500
ROOT_DIR="KEGG/500_size"
runOther PDX2_TET1_g1g2g3.sig.down.txt $ROOT_DIR ENSEMBL "PDX2_TET1_g1g2g3_sig_down" 500

RESULT_DIR="GO/500_size"
runGO PDX2_TET1_g1g2g3.sig.up.txt $RESULT_DIR ENSEMBL "PDX2_TET1_g1g2g3_sig_up" 500
ROOT_DIR="KEGG/500_size"
runOther PDX2_TET1_g1g2g3.sig.up.txt $ROOT_DIR ENSEMBL "PDX2_TET1_g1g2g3_sig_up" 500

RESULT_DIR="GO/500_size"
runGO PDX2_TET1_g1g2g3.sig.txt $RESULT_DIR ENSEMBL "PDX2_TET1_g1g2g3_sig" 500
ROOT_DIR="$ENRICH_DIR"
runOther PDX2_TET1_g1g2g3.sig.txt $ROOT_DIR ENSEMBL "PDX2_TET1_g1g2g3_sig" 500
