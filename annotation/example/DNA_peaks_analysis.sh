#/usr/bin/sh
BASE=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/hyweng/tet1/chip_seq/publicData

# functions
function AnnoBed {
  BED_INPUT_ARR=$1
  BED_OUTPUT_ARR=$2
  SIGNAL_VAL=$3
  PVAL=$4
  for i in "${!BED_INPUT_ARR[@]}";
  do
    originalBed=${BED_INPUT_ARR[$i]}
    bed=${BED_OUTPUT_ARR[$i]}
    awk -v signal="$SIGNAL_VAL" -v pval="$PVAL" \
      'BEGIN{FS="\t";OFS="\t";}
      {
        if(/^#/){
          print
        }else{
          if ($7>=signal && $8>=pval){
            print
          }
        }
      }' $originalBed > $bed
    baseName=${bed%.*}
    peakAnnotate.py -input $bed \
      --keepName -mode DNA -gsize $GENOME_SIZE \
      -anno $FULL_BED -geneClassFile $GENETYPE_FILE \
      -extraAnno $PRIMARY_MIRNA_BED $T_RNA $R_RNA $CPG $STATELLINE $REPEAT \
      -extraType gene gene gene element element element\
      -ncbiGeneInfo $NCBI_GENE_INFO \
      -output ${baseName}.anno.txt
  done
}

function FeatureStats {
  annoResFold=$1
  statsFile=$2
  peakType=$3
  echo -e "ExpName\tExpType\tPeakType\tFeature\tCount\tPercentage" > $statsFile
  for annoResFile in `find $annoResFold -type f -name "*.anno.txt" | sort`;
  do
    baseName=${annoResFile%%.anno.txt}
    baseName=${baseName##*/}
    expName=${baseName%%_control*}
    if [[ $annoResFile =~ 'cutoff' ]]; then
      expType="cutoff";
    else
      expType="original";
    fi
    awk -v expName="$expName" -v expType="$expType" -v peakType="$peakType" \
      'BEGIN{FS="\t";OFS="\t";sum=0;}
      {
        if (FNR >1 ){
          if (peakType == "narrowPeak") {
            mainFeature = $21;
            extraName = $27;
            extraFeature = $28;
          }else{
            mainFeature = $20;
            extraName = $26;
            extraFeature = $27;
          }
          if (mainFeature == "intergenic") {
            if (extraFeature == "na" || extraFeature == "intergenic") {
              feature="intergenic";
            }else {
              feature=extraFeature;
              if (extraFeature == "repeat" && extraName == "Alu"){
                feature="Alu";
              }
            }
          }else{
            feature=mainFeature;
          }
          featureArr[feature]+=1
          sum += 1
        }
      }
      END{
        for(i=1;i<=asorti(featureArr,key);i++){
          feature = key[i];
          fCount = featureArr[feature];
          percentage = sprintf("%.2f", fCount/sum*100);
          print expName, expType, peakType, feature, fCount, percentage;
        }
      }' $annoResFile >> $statsFile
  done
}

function GeneTypeStats {
  annoResFold=$1
  statsFile=$2
  peakType=$3
  echo -e "ExpName\tExpType\tPeakType\tGeneType\tCount\tPercentage" > $statsFile
  for annoResFile in `find $annoResFold -type f -name "*.anno.txt" | sort`;
  do
    baseName=${annoResFile%%.anno.txt}
    baseName=${baseName##*/}
    expName=${baseName%%_control*}
    if [[ $annoResFile =~ 'cutoff' ]]; then
      expType="cutoff";
    else
      expType="original";
    fi
    awk -v expName="$expName" -v expType="$expType" -v peakType="$peakType" \
      'BEGIN{FS="\t";OFS="\t";sum=0;}
      {
        if (FNR >1 ){
          if (peakType == "narrowPeak") {
            mainGeneType = $16;
            mainGeneClass = $17;
            extraName = $27;
            extraFeature = $28;
            extraGeneType = $29;
          }else{
            mainGeneType = $15;
            mainGeneClass = $16;
            extraName = $26;
            extraFeature = $27;
            extraGeneType = $28;
          }
          if (mainGeneType == "protein_coding") {
            mainGeneClass = "protein_coding";
          }
          if (mainGeneClass == "intergenic") {
            if (extraFeature == "na" || extraFeature == "intergenic") {
              geneType="intergenic";
            }else {
              geneType=extraGeneType;
              if (extraFeature == "repeat" && extraName == "Alu"){
                geneType="Alu";
              }
            }
          }else{
            geneType=mainGeneClass;
          }
          geneTypeArr[geneType]+=1
          sum += 1;
        }
      }
      END{
        for(i=1;i<=asorti(geneTypeArr,key);i++){
          feature = key[i];
          fCount = geneTypeArr[feature];
          percentage = sprintf("%.2f", fCount/sum*100);
          print expName, expType, peakType, feature, fCount, percentage;
        }
      }' $annoResFile >> $statsFile
  done
}


# analysis start
## for mouse
## genome
PUBLIC_BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome"
GENOME_SIZE="$PUBLIC_BASE/size/mm10.chrom.sizes"
GENETYPE_FILE="$PUBLIC_BASE/annotation/mm10/geneType.mm10.txt"
## annotation
NCBI_GENE_INFO="$PUBLIC_BASE/annotation/mm10/mm10.ncbi.gene_info"
FULL_BED="$PUBLIC_BASE/annotation/mm10/vM24/gencode.vM24.annotation.anno.bed12"
PRIMARY_MIRNA_BED="$PUBLIC_BASE/annotation/mm10/addition/miRBase.v22.primary.bed"
CPG="$PUBLIC_BASE/annotation/mm10/addition/cpgisland.bed"
STATELLINE="$PUBLIC_BASE/annotation/mm10/addition/microsatelline.bed"
REPEAT="$PUBLIC_BASE/annotation/mm10/addition/repeat.bed"
R_RNA="$PUBLIC_BASE/annotation/mm10/addition/rRNA.bed"
T_RNA="$PUBLIC_BASE/annotation/mm10/addition/tRNA.bed"
########### mapping & calling ##########
########### annotation ##########

ANALYSIS_DIR="$BASE/analysis/mouse"
cd $ANALYSIS_DIR

## delete previous results
find $ANALYSIS_DIR -type f | grep -v "original" | xargs -I {} rm -f {}

# original peaks
for fold in `find $ANALYSIS_DIR -type d | grep "original"`;
do
  PEAK_DIR=`realpath $fold`
  ANNO_RESULT_DIR="${PEAK_DIR}/../"
  declare -a BED_INPUT_ARR=()
  declare -a BED_OUTPUT_ARR=()
  for peak in `find $fold -type f`;
  do
    peakFileName=$(basename $peak)
    extension="${peakFileName##*.}"
    baseName=$(basename $peakFileName ".$extension")
    BED_INPUT_ARR+=( "$peak" )
    BED_OUTPUT_ARR+=( "${ANNO_RESULT_DIR}/${peakFileName}" )
  done
  AnnoBed $BED_INPUT_ARR $BED_OUTPUT_ARR 0 0
done

# cutoff filter with fold_enrichmen>=2, p<=1e-4
for fold in `find $ANALYSIS_DIR -type d | grep "original"`;
do
  PEAK_DIR=`realpath $fold`
  ANNO_RESULT_DIR="${PEAK_DIR}/../"
  declare -a BED_INPUT_ARR=()
  declare -a BED_OUTPUT_ARR=()
  for peak in `find $fold -type f`;
  do
    peakFileName=$(basename $peak)
    extension="${peakFileName##*.}"
    baseName=$(basename $peakFileName ".$extension")
    BED_INPUT_ARR+=( "$peak" )
    BED_OUTPUT_ARR+=( "${ANNO_RESULT_DIR}/${baseName}.cutoff.${extension}" )
  done
  AnnoBed $BED_INPUT_ARR $BED_OUTPUT_ARR 2 4
done

###
# feature statistics
FeatureStats $ANALYSIS_DIR/broadPeak broadPeakFeatureStats.stats broadPeak
FeatureStats $ANALYSIS_DIR/narrowPeak narrowPeakFeatureStats.stats narrowPeak

GeneTypeStats $ANALYSIS_DIR/broadPeak broadPeakGeneTypeStats.stats broadPeak
GeneTypeStats $ANALYSIS_DIR/narrowPeak narrowPeakGeneTypeStats.stats narrowPeak
