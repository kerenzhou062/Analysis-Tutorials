#!/bin/sh

BASE="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/hyweng/tet1/chip_seq/publicData"
PEAK_DIR="$BASE/peak"
ANALYSIS_DIR="$BASE/analysis/mouse"
TET1_PEAK_DIR="/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/chshen/tet1/publicData/peak"
## genome

if [[ -d $ANALYSIS_DIR ]]; then
  mkdir -p $ANALYSIS_DIR
fi

cd $ANALYSIS_DIR

# make directories
narrowPeakDir="$ANALYSIS_DIR/narrowPeak/original"
broadPeakDir="$ANALYSIS_DIR/broadPeak/original"
rm -rf $narrowPeakDir
rm -rf $broadPeakDir
mkdir -p $narrowPeakDir
mkdir -p $broadPeakDir

# copy peak files
find $PEAK_DIR -type f -name "*.narrowPeak" | grep -P -v 'pr1|pr2|with\w*?Pr|IDR' | grep "mouse" | \
  xargs -I {} cp {} $narrowPeakDir

find $PEAK_DIR -type f -name "*.broadPeak" | grep -P -v 'pr1|pr2|with\w*?Pr|IDR' | grep "mouse" | \
  xargs -I {} cp {} $broadPeakDir

find ./ -type f | grep "HSC_H4K16ac" | grep "_pooled_peaks" | xargs -I {} rm -f {}

find $TET1_PEAK_DIR -type f -name "*.narrowPeak" | grep -P -v 'pr1|pr2|with\w*?Pr|IDR' | grep -v "_pooled_peaks" | \
  xargs -I {} cp {} $narrowPeakDir

find $TET1_PEAK_DIR -type f -name "*.broadPeak" | grep -P -v 'pr1|pr2|with\w*?Pr|IDR' | grep -v "_pooled_peaks" | \
  xargs -I {} cp {} $broadPeakDir

# intersect 
for i in $narrowPeakDir/E14_H4K16ac*;
do
  extension="${i##*.}"
  baseName=$(basename $i ".$extension")
  bedtools intersect -a $narrowPeakDir/mES_Tet1_WT_rep1_final_peaks.${extension} \
    -b $i > $narrowPeakDir/tet1_intersect_${baseName}.${extension}
done

# intersect 
for i in $broadPeakDir/E14_H4K16ac*;
do
  extension="${i##*.}"
  baseName=$(basename $i ".$extension")
  bedtools intersect -a $broadPeakDir/mES_Tet1_WT_rep1_final_peaks.${extension} \
    -b $i > $broadPeakDir/tet1_intersect_${baseName}.${extension}
done
#

for i in `find ./ -type f -name "*.narrowPeak"`;
do
  sort -t $'\t' -k 1,1 -k2,2n -o $i $i
  sed -i '1 i\#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tenrichment\tpValue(-log10)\tqValue(-log10)\tpeak' $i
done

for i in `find ./ -type f -name "*.broadPeak"`;
do
  sort -t $'\t' -k 1,1 -k2,2n -o $i $i
  sed -i '1 i\#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tenrichment\tpValue(-log10)\tqValue(-log10)' $i
done

