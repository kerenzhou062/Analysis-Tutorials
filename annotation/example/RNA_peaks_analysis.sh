#/usr/bin/sh
BASE=/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/brandontan/ythdc1/public_data
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
R_RNA="$PUBLIC_BASE/annotation/hg38/addition/rRNA.bed"
T_RNA="$PUBLIC_BASE/annotation/hg38/addition/tRNA.bed"
REPEAT="$PUBLIC_BASE/annotation/hg38/addition/repeat.bed"

########### annotation ##########

function annoBed {
  BED_INPUT_ARR=$1
  BED_OUTPUT_ARR=$2
  BED_NAME_ARR=$3
  for i in "${!BED_INPUT_ARR[@]}";
  do
    originalBed=$FINAL_DIR/${BED_INPUT_ARR[$i]}
    bed=$CUSTOM_DIR/${BED_OUTPUT_ARR[$i]}
    awk -v p="${BED_NAME_ARR[$i]}" \
      'BEGIN{FS="\t";OFS="\t";}
      {
        if(/^#/){
          print
        }else{
          $4=p""FNR;print
        }
      }' $originalBed > $bed
    baseName=${bed%%.bed}
    peakAnnotate.py -input $bed \
      -strand -name "${BED_NAME_ARR[$i]}" -mode RNA -codon "100,100" \
      -anno $FULL_BED -geneClassFile $GENETYPE_FILE \
      -extraAnno $SNORNA_BED $PRIMARY_MIRNA_BED $R_RNA $T_RNA $REPEAT \
      -extraType gene gene gene gene element\
      -ncbiGeneInfo $NCBI_GENE_INFO \
      -output ${baseName}.anno.txt
  done
}


function peakFeatureStats {
  annoResFile=$1
  expName=$2
  expType=$3
  statsFile=$4
  awk -v name="$expName" -v type="$expType" \
    'BEGIN{FS="\t";OFS="\t";sum=0;}
    {
      if (FNR >1 ){
        if ($(NF-10) == "intergenic") {
          if ($(NF-2) == "na") {
            feature="intergenic";
          }else {
            feature=$(NF-2);
            if ($(NF-2) == "repeat" && $(NF-3) == "Alu"){
              feature="Alu";
            }
          }
        }else{
          feature=$(NF-10);
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
    }' $annoResFile >> $statsFile
}

function peakGeneTypeStats {
  annoResFile=$1
  expName=$2
  expType=$3
  statsFile=$4
  awk -v name="$expName" -v type="$expType" \
    'BEGIN{FS="\t";OFS="\t";sum=0;}
    {
      if (FNR >1 ){
        if ($(NF-14) == "intergenic") {
          if ($(NF-1) == "na") {
            geneType="intergenic";
          }else {
            geneType=$(NF-1);
            if ($(NF-1) == "repeat" && $(NF-3) == "Alu"){
              geneType="Alu";
            }
          }
        }else{
          geneType=$(NF-14);
          if ($(NF-1) == "snoRNA") {
            geneType=$(NF-1);
          }else if ($(NF-1) == "miRNA_primary_transcript") {
            geneType=$(NF-1);
          }
        }
        geneTypeArr[geneType]+=1
      }
      sum += 1
    }
    END{
      for(i=1;i<=asorti(geneTypeArr,key);i++){
        geneType = key[i];
        fCount = geneTypeArr[geneType];
        percentage = sprintf("%.2f", fCount/sum*100);
        print name, type, geneType, fCount, percentage;
      }
    }' $annoResFile >> $statsFile
}

## CTK_BWA
if [[ ! -d $CUSTOM_DIR ]]; then
  mkdir -p $CUSTOM_DIR
fi

cd $CUSTOM_DIR

## CTK_BWA
FINAL_DIR=$BASE/CTK_BWA/final
BED_INPUT_ARR=( "HeLa_YTHDC1_WT-1_IP.peak.sig.t2c.bed" "HeLa_YTHDC1_WT-1_IP.peak.sig.t2c.RRACH.bed" )
BED_INPUT_ARR+=( "HeLa_YTHDC1_WT-1_IP.t2c.CIMS.RRACH.bed" )
BED_OUTPUT_ARR=( "HeLa_YTHDC1_WT-1_IP.CTK.peak.sig.t2c.bed" "HeLa_YTHDC1_WT-1_IP.CTK.peak.sig.t2c.RRACH.bed" )
BED_OUTPUT_ARR+=( "HeLa_YTHDC1_WT-1_IP.CTK.t2c.CIMS.RRACH.bed" )
BED_NAME_ARR=( "YTHDC1|CTK_BWA|peak|sig|t2c=" "YTHDC1|CTK_BWA|peak|sig|t2c|RRACH=" "YTHDC1|CTK_BWA|CIMS|RRACH=" )
#annoBed $BED_INPUT_ARR $BED_OUTPUT_ARR $BED_NAME_ARR

## PARalyzer_bowtie
FINAL_DIR=$BASE/PARalyzer_bowtie/final
BED_INPUT_ARR=( "HeLa_YTHDC1_WT-1_IP.cluster.bed" "HeLa_YTHDC1_WT-1_IP.cluster.RRACH.bed" )
BED_OUTPUT_ARR=( "HeLa_YTHDC1_WT-1_IP.PARalyzer.bowtie.cluster.bed" "HeLa_YTHDC1_WT-1_IP.PARalyzer.bowtie.cluster.RRACH.bed" )
BED_NAME_ARR=( "YTHDC1|PARalyzer_bowtie|cluster=" "YTHDC1|PARalyzer_bowtie|cluster|RRACH=" )
#annoBed $BED_INPUT_ARR $BED_OUTPUT_ARR $BED_NAME_ARR

## PARalyzer_STAR
FINAL_DIR=$BASE/PARalyzer_STAR/final
BED_INPUT_ARR=( "HeLa_YTHDC1_WT-1_IP.cluster.bed" "HeLa_YTHDC1_WT-1_IP.cluster.RRACH.bed" )
BED_OUTPUT_ARR=( "HeLa_YTHDC1_WT-1_IP.PARalyzer.STAR.cluster.bed" "HeLa_YTHDC1_WT-1_IP.PARalyzer.STAR.cluster.RRACH.bed" )
BED_NAME_ARR=( "YTHDC1|PARalyzer_STAR|cluster=" "YTHDC1|PARalyzer_STAR|cluster|RRACH=" )
#annoBed $BED_INPUT_ARR $BED_OUTPUT_ARR $BED_NAME_ARR

# CTK
bedtools intersect -a HeLa_YTHDC1_WT-1_IP.CTK.peak.sig.t2c.bed -b HeLa_YTHDC1_WT-1_IP.PARalyzer.bowtie.cluster.RRACH.bed -s -u | \
  bedtools intersect -a stdin -b HeLa_YTHDC1_WT-1_IP.PARalyzer.STAR.cluster.RRACH.bed -s -u | wc -l >
echo -e "ExpName\tExpType\tFeature\tCount\tPercentage" > peakpeakFeatureStats.stats
echo -e "ExpName\tExpType\tGeneTyp\tCount\tPercentage" > peakpeakGeneTypeStats.stats
for annoResFile in `find $CUSTOM_DIR -maxdepth 1 -type f -name "*.anno.txt" | grep -v "RRACH" |sort`;
do
  baseName=${annoResFile%%.anno.txt}
  baseName=${baseName##*/}
  expName=${baseName%%_IP*}
  if [[ $annoResFile =~ 'PARalyzer.bowtie' ]]; then
    expType="PARalyzer_bowtie"
  elif  [[ $annoResFile =~ 'CTK' ]]; then
    expType="CTK_BWA"
  elif  [[ $annoResFile =~ 'PARalyzer.STAR' ]]; then
    expType="PARalyzer_STAR"
  fi
  peakFeatureStats $annoResFile $expName $expType peakpeakFeatureStats.stats
  peakGeneTypeStats $annoResFile $expName $expType peakpeakGeneTypeStats.stats
done

