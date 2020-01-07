#/usr/bin/sh
PIPELINE_DIR=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/home/public/pipeline/PAR-CLIP

BASE=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/hyweng/tet1
CUSTOM_DIR="$BASE/custom_analysis"
## genome
PUBLIC_BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/public/genome"
GENOME_SIZE="$PUBLIC_BASE/size/number/hg38.chrom.sizes"
GENETYPE_FILE="$PUBLIC_BASE/annotation/hg38/geneType.hg38.txt"
## annotation
NCBI_GENE_INFO="$PUBLIC_BASE/annotation/hg38/hg38.ncbi.gene_info"
FULL_BED="$PUBLIC_BASE/annotation/hg38/v32/gencode.v32.annotation.anno.bed12"
SNORNA_BED="$PUBLIC_BASE/annotation/hg38/addition/snoRNA.bed"
PRIMARY_MIRNA_BED="$PUBLIC_BASE/annotation/hg38/addition/miRBase.v22.primary.bed"
CPG="$PUBLIC_BASE/annotation/hg38/addition/cpgisland.bed"
STATELLINE="$PUBLIC_BASE/annotation/hg38/addition/microsatelline.bed"
REPEAT="$PUBLIC_BASE/annotation/hg38/addition/repeat.bed"
R_RNA="$PUBLIC_BASE/annotation/hg38/addition/rRNA.bed"
T_RNA="$PUBLIC_BASE/annotation/hg38/addition/tRNA.bed"

########### function ##########

function annoBed {
  BED_INPUT_ARR=$1
  BED_OUTPUT_ARR=$2
  BED_NAME_ARR=$3
  CUTOFF=$4
  for i in "${!BED_INPUT_ARR[@]}";
  do
    originalBed=$PEAK_DIR/${BED_INPUT_ARR[$i]}
    bed=$CUSTOM_DIR/${BED_OUTPUT_ARR[$i]}
    awk -v prefix="${BED_NAME_ARR[$i]}" -v cutoff="$CUTOFF" \
      'BEGIN{FS="\t";OFS="\t";}
      {
        if(/^#/){
          print
        }else{
          if(cutoff==1){
            if ($8>=4 && $7>=2){
              $4=prefix""FNR;print
            }
          }else{
            $4=prefix""FNR;print
          }
        }
      }' $originalBed > $bed
    baseName=${bed%%.narrowPeak}
    peakAnnotate.py -input $bed \
      -name "Kas1|macs2=" -mode DNA -gsize $GENOME_SIZE \
      -anno $FULL_BED -geneClassFile $GENETYPE_FILE \
      -extraAnno $SNORNA_BED $PRIMARY_MIRNA_BED $T_RNA $R_RNA $CPG $STATELLINE $REPEAT \
      -extraType gene gene gene gene element element element\
      -ncbiGeneInfo $NCBI_GENE_INFO \
      -output ${baseName}.anno.txt
  done
}

function peakFeatureStats {
  expName=$1
  expType=$2
  input=$3
  output=$4
  awk -v name="$expName" -v type="$expType" \
    'BEGIN{FS="\t";OFS="\t";sum=0;}
    {
      if (FNR >1 ){
        if ($21 == "intergenic") {
          if ($28 == "na" || $28 == "intergenic") {
            feature="intergenic";
          }else {
            feature=$28;
            if ($28 == "repeat" && $27 == "Alu"){
              feature="Alu";
            }
          }
        }else{
          feature=$21;
        }
        featureArr[feature]+=1
      }
      sum += 1
    }
    END{
      for(i=1;i<=asorti(featureArr,key);i++){
        feature = key[i];
        fCount = featureArr[feature];
        percentage = sprintf("%.2f", fCount/sum*100);
        print name, type, feature, fCount, percentage;
      }
    }' $input >> $output
}

function peakGeneTypeStats {
  expName=$1
  expType=$2
  input=$3
  output=$4
  awk -v name="$expName" -v type="$expType" \
    'BEGIN{FS="\t";OFS="\t";sum=0;}
    {
      if (FNR >1 ){
        if ($17 == "intergenic") {
          if ($28 == "na" || $28 == "intergenic") {
            geneType="intergenic";
          }else {
            geneType=$29;
            if ($28 == "repeat" && $27 == "Alu"){
              geneType="Alu";
            }
          }
        }else{
          geneType=$17;
        }
        geneTypeArr[geneType]+=1
      }
      sum += 1;
    }
    END{
      for(i=1;i<=asorti(geneTypeArr,key);i++){
        geneType = key[i];
        fCount = geneTypeArr[geneType];
        percentage = sprintf("%.2f", fCount/sum*100);
        print name, type, geneType, fCount, percentage;
      }
    }' $input >> $output
}

########### annotation ##########

if [[ ! -d $CUSTOM_DIR ]]; then
  mkdir -p $CUSTOM_DIR
fi

cd $CUSTOM_DIR

# original
PEAK_DIR=$CUSTOM_DIR/original
BED_INPUT_ARR=( "Kas1_control_NA_rep1_final_peaks.narrowPeak" "MA93_control_NA_rep1_final_peaks.narrowPeak")
BED_INPUT_ARR+=( "MDSL_control_NA_rep1_final_peaks.narrowPeak" )
BED_OUTPUT_ARR=( "Kas1_control_NA_rep1_final_peaks.narrowPeak" "MA93_control_NA_rep1_final_peaks.narrowPeak")
BED_OUTPUT_ARR+=( "MDSL_control_NA_rep1_final_peaks.narrowPeak" )
BED_NAME_ARR=( "Kas1|macs2|peak=" "MA93|macs2|peak=" "MDSL|macs2|peak=" )
annoBed $BED_INPUT_ARR $BED_OUTPUT_ARR $BED_NAME_ARR 0

# filter with p<1e-4, fold_enrichmen>=2
PEAK_DIR=$CUSTOM_DIR/original
BED_INPUT_ARR=( "Kas1_control_NA_rep1_final_peaks.narrowPeak" "MA93_control_NA_rep1_final_peaks.narrowPeak")
BED_INPUT_ARR+=( "MDSL_control_NA_rep1_final_peaks.narrowPeak" )
BED_OUTPUT_ARR=( "Kas1_control_NA_rep1_final_peaks.cutoff.narrowPeak" "MA93_control_NA_rep1_final_peaks.cutoff.narrowPeak")
BED_OUTPUT_ARR+=( "MDSL_control_NA_rep1_final_peaks.cutoff.narrowPeak" )
BED_NAME_ARR=( "Kas1|macs2|peak|cutoff=" "MA93|macs2|peak|cutoff=" "MDSL|macs2|peak|cutoff=" )
annoBed $BED_INPUT_ARR $BED_OUTPUT_ARR $BED_NAME_ARR 1

###
# feature statistics & geneType statistics
echo -e "ExpName\tExpType\tFeature\tCount\tPercentage" > peakFeatureStats.stats
echo -e "ExpName\tExpType\tGeneType\tCount\tPercentage" > peakpeakGeneTypeStats.stats

for annoResFile in `find $CUSTOM_DIR -maxdepth 1 -type f -name "*.anno.txt" | sort`;
do
  baseName=${annoResFile%%.anno.txt}
  baseName=${baseName##*/}
  expName=${baseName%%_control*}
  if [[ $annoResFile =~ 'cutoff' ]]; then
    expType="cutoff";
  else
    expType="original";
  fi
  peakFeatureStats $expName $expType $annoResFile peakFeatureStats.stats
  peakGeneTypeStats $expName $expType $annoResFile peakpeakGeneTypeStats.stats
done
