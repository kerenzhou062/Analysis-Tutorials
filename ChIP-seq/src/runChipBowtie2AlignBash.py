#!/usr/bin/env python3
import os
import sys
import argparse
import re
from glob import glob
from collections import defaultdict
import datetime
from PubAlbum import Anno

#usage: runChipSeBowtie2AlignBash.py or runChipSeBowtie2AlignBash.py <fastq dir>

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-blacklist', action='store', type=str,
                    default='hg38', help='blacklist for genome: prefix or file (eg. hg38)')
parser.add_argument('-cpu', action='store', type=str,
                    default='10', help='threads used for bowtie2')
parser.add_argument('-fasta', action='store', type=str,
                    default='hg38', help='fasta for genome: prefix or file (eg. hg38)')
parser.add_argument('-input', action='store', type=str,
                    default='./', help='input fastq directory (eg. HepG2_shWTAP_IP_rep1_run1_2.fastq)')
parser.add_argument('-index', action='store', type=str,
                    default='hg38', help='bowtie2 genome index: prefix or file (eg. hg38)')
parser.add_argument('-memory', action='store', type=str,
                    default='50G', help='memory used for sbatch')
parser.add_argument('-quality', action='store', type=str,
                    default='30', help='quality for filtering fastq')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    print("Running runChipSeBowtie2AlignBash.py with defaultdict parameters...")

#public arguments
threadNum = args.cpu
memory = args.memory
quality = args.cpu
anno = Anno()
genomeIndex = anno.gindex(args.index)
blacklist = anno.blacklist(args.blacklist)
fasta = anno.fasta(args.fasta)
basepath = os.path.realpath(args.input)

#filename: HepG2_shWTAP_IP_rep1_run1_1.fastq, HepG2_shWTAP_IP_rep1_run1_2.fastq, HepG2_shWTAP_IP_run_rep1_NA.fastq
extsList = ["*.fastq", '*.fq', '*.fastq.gz', '*.fq.gz']
fastqList = sorted([f for ext in extsList for f in glob(os.path.join(basepath, ext))])
bashDir = os.path.join(basepath, 'runBash')
logDir = os.path.join(bashDir, 'log')
os.makedirs(logDir, exist_ok=True)

readLen = 0
fqExt = '' # .fq.gz|.fq|.fastq|.fastq.gz
expFqDict = defaultdict(dict)
# expFqDict: exp->IP->rep->[fq1,fq2]
for fastq in fastqList:
    basename = os.path.basename(fastq)
    fqExt = re.findall(r'\..+$', basename)[0]
    cleanBasename = re.sub(r'\..+$', '', basename)
    #HepG2 shWTAP IP rep1 1
    basenameList = cleanBasename.split('_')
    #HepG2_shWTAP
    exp = '_'.join(basenameList[0:2])
    ip = basenameList[2]
    rep = basenameList[3]
    run = basenameList[4]
    pairedNum = basenameList[5]
    if ip in expFqDict[exp]:
        if pairedNum in expFqDict[exp][ip][rep]:
            expFqDict[exp][ip][rep][pairedNum].append(basename)
        else:
            expFqDict[exp][ip][rep][pairedNum] = [basename]
    else:
        expFqDict[exp][ip] = defaultdict(dict)
        expFqDict[exp][ip][rep][pairedNum] = [basename]
    if readLen == 0:
        if re.search(r'gz', fqExt):
            import gzip
            with gzip.open(fastq,'rt') as f:
                for i, line in enumerate(f):
                    if i == 1:
                        readLen = len(line.strip())
                    elif i > 1:
                        break
        else:
            with open(fastq,'rt') as f:
                for i, line in enumerate(f):
                    if i == 1:
                        readLen = len(line.strip())
                    elif i > 1:
                        break

trim = 50
if readLen < trim:
    trim = readLen

mainAlignDir = os.path.join(basepath, 'alignment')

# ================================
#step template
# ================================
#mapping command
step1aSeCommand = '''
# mapping single-end reads
zcat -f ${{fastqR1}} | bowtie2 --mm -x ${{bwt2_idx}} --threads {threadNum} -U - 2> ${{log}} | \\
  samtools view -Su /dev/stdin | \\
  samtools sort -T ${{prefix}} -O bam -o ${{prefix}}.bam

'''

setp1aPeCommand = '''
# mapping paired-end reads
bowtie2 -X2000 --mm -x ${{bwt2_idx}} --threads {threadNum} \
  -1 $fastqR1 -2 $fastqR2 2> ${{log}} | \\
  samtools view -Su /dev/stdin | \\
  samtools sort -T ${{prefix}} -O bam -o ${{prefix}}.bam

'''

#post-alignment filtering 2a steps
step1bSeCommand = '''
OFPREFIX="${{prefix}}"
RAW_BAM_FILE="${{OFPREFIX}}.bam"
FILT_BAM_PREFIX="${{OFPREFIX}}.filt.srt"
FILT_BAM_FILE="${{FILT_BAM_PREFIX}}.bam"
MAPQ_THRESH={quality}

samtools view -F 1804 -q ${{MAPQ_THRESH}} -b ${{RAW_BAM_FILE}} -o ${{FILT_BAM_FILE}}
# ======================================
# Mark duplicates
# ====================================
echo "Mark duplicates..."

TMP_FILT_BAM_FILE="${{FILT_BAM_PREFIX}}.dupmark.bam"
MARKDUP="/opt/picard/2.21.1/picard.jar MarkDuplicates"
DUP_FILE_QC="${{FILT_BAM_PREFIX}}.dup.qc" # QC file

java -Xmx4G -jar ${{MARKDUP}} INPUT=${{FILT_BAM_FILE}} \\
  OUTPUT=${{TMP_FILT_BAM_FILE}} \\
  METRICS_FILE=${{DUP_FILE_QC}} \\
  VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

mv ${{TMP_FILT_BAM_FILE}} ${{FILT_BAM_FILE}}

# ==========================================
# Remove duplicates
# Index final position sorted BAM
# ==========================================
echo "Remove duplicates and index BAM..."

FINAL_BAM_PREFIX="${{OFPREFIX}}.filt.srt.nodup"
FINAL_BAM_FILE="${{FINAL_BAM_PREFIX}}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${{FINAL_BAM_PREFIX}}.bai" # To be stored
FINAL_BAM_FILE_MAPSTATS="${{FINAL_BAM_PREFIX}}.flagstat.qc" # QC file

# Index Final BAM file
samtools view -F 1804 -b ${{FILT_BAM_FILE}} -o ${{FINAL_BAM_FILE}}
samtools index ${{FINAL_BAM_FILE}} ${{FINAL_BAM_INDEX_FILE}}
samtools sort -n --threads 10 ${{FINAL_BAM_FILE}} -O SAM | \
  SAMstats --sorted_sam_file - --outf ${{FINAL_BAM_FILE_MAPSTATS}}

# ================================
# Compute library complexity
# ================================
# sort by position and strand
# Obtain unique count statistics
echo "Compute library complexity..."

# PBC File output
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab]
# PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
PBC_FILE_QC="${{FINAL_BAM_PREFIX}}.pbc.qc"

bedtools bamtobed -i ${{FILT_BAM_FILE}} | \\
  awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$6}}' | \\
  grep -v 'chrM' | sort | uniq -c | \\
  awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}}
    END{{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}}' \\
  > ${{PBC_FILE_QC}}

# do not delete ${{FILT_BAM_FILE}}
## rm ${{FILT_BAM_FILE}}

# link ${{FINAL_BAM_FILE}} to ${{ipDir}}
ln -sf ${{FINAL_BAM_FILE}} ${{ipDir}}

'''

step1bPeCommand = '''
OFPREFIX="${{prefix}}"
RAW_BAM_FILE="${{OFPREFIX}}.bam"
FILT_BAM_PREFIX="${{OFPREFIX}}.filt.srt"
FILT_BAM_FILE="${{FILT_BAM_PREFIX}}.bam"
TMP_FILT_BAM_PREFIX="${{FILT_BAM_PREFIX}}.tmp.nmsrt"
TMP_FILT_BAM_FILE="${{TMP_FILT_BAM_PREFIX}}.bam"
MAPQ_THRESH={quality}

# Will produce name sorted BAM
samtools view -F 1804 -f 2 -q ${{MAPQ_THRESH}} -u ${{RAW_BAM_FILE}} | \\
  samtools sort -n - -T ${{TMP_FILT_BAM_PREFIX}} -O bam -o ${{TMP_FILT_BAM_FILE}}

# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM
samtools fixmate -r ${{TMP_FILT_BAM_FILE}} ${{OFPREFIX}}.fixmate.tmp
samtools view -F 1804 -f 2 -u ${{OFPREFIX}}.fixmate.tmp | \\
  samtools sort - -T ${{FILT_BAM_PREFIX}} -O bam -o ${{FILT_BAM_FILE}}
rm ${{OFPREFIX}}.fixmate.tmp
rm ${{TMP_FILT_BAM_FILE}}

# ======================================
# Mark duplicates
# ====================================
echo "Mark duplicates..."

TMP_FILT_BAM_FILE="${{FILT_BAM_PREFIX}}.dupmark.bam"
MARKDUP="/opt/picard/2.21.1/picard.jar MarkDuplicates"
DUP_FILE_QC="${{FILT_BAM_PREFIX}}.dup.qc"
java -Xmx4G -jar ${{MARKDUP}} INPUT=${{FILT_BAM_FILE}} \\
  OUTPUT=${{TMP_FILT_BAM_FILE}} \\
  METRICS_FILE=${{DUP_FILE_QC}} \\
  VALIDATION_STRINGENCY=LENIENT \\
  ASSUME_SORTED=true REMOVE_DUPLICATES=false

mv ${{TMP_FILT_BAM_FILE}} ${{FILT_BAM_FILE}}

# ==========================================
# Remove duplicates
# Index final position sorted BAM
# Create final name sorted BAM
# ==========================================
echo "Remove duplicates and index BAM..."

FINAL_BAM_PREFIX="${{OFPREFIX}}.filt.srt.nodup"
FINAL_BAM_FILE="${{FINAL_BAM_PREFIX}}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${{FINAL_BAM_PREFIX}}.bai"
FINAL_BAM_FILE_MAPSTATS="${{FINAL_BAM_PREFIX}}.flagstat.qc" # QC file
FINAL_NMSRT_BAM_PREFIX="${{OFPREFIX}}.filt.nmsrt.nodup"
FINAL_NMSRT_BAM_FILE="${{FINAL_NMSRT_BAM_PREFIX}}.bam" # To be stored

samtools view -F 1804 -f 2 -b ${{FILT_BAM_FILE}} -o ${{FINAL_BAM_FILE}}
samtools sort -n ${{FINAL_BAM_FILE}} -T ${{FINAL_NMSRT_BAM_PREFIX}} -o ${{FINAL_NMSRT_BAM_FILE}}
# Index Final BAM file
samtools index ${{FINAL_BAM_FILE}} ${{FINAL_BAM_INDEX_FILE}}
samtools sort -n --threads 10 ${{FINAL_BAM_FILE}} -O SAM | \\
  SAMstats --sorted_sam_file - --outf ${{FINAL_BAM_FILE_MAPSTATS}}

# ================================
# Compute library complexity
# ================================
# Sort by name
# convert to bedPE and obtain fragment coordinates
# sort by position and strand
# Obtain unique count statistics
echo "Compute library complexity..."

PBC_FILE_QC="${{FINAL_BAM_PREFIX}}.pbc.qc"
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab]
# NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

samtools sort -n ${{FILT_BAM_FILE}} -o ${{OFPREFIX}}.srt.tmp.bam
bedtools bamtobed -bedpe -i ${{OFPREFIX}}.srt.tmp.bam | \\
  awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$4,$6,$9,$10}}' | grep -v 'chrM' | sort | uniq -c | \\
  awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}}
    ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} 
    END{{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}}' \\
  > ${{PBC_FILE_QC}}
rm ${{OFPREFIX}}.srt.tmp.bam

# do not delete ${{FILT_BAM_FILE}}
## rm ${{FILT_BAM_FILE}}

# link ${{FINAL_BAM_FILE}} to ${{ipDir}}
ln -sf ${{FINAL_BAM_FILE}} ${{ipDir}}

'''

step2aSeCommand = '''
FINAL_TA_FILE="${{FINAL_BAM_PREFIX}}.final.tagAlign.gz"

bedtools bamtobed -i ${{FINAL_BAM_FILE}} | \
  awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}' | \\
  gzip -nc > ${{FINAL_TA_FILE}}

# link ${{FINAL_TA_FILE}} to ${{ipDir}}
ln -sf ${{FINAL_TA_FILE}} ${{ipDir}}

'''

step2aPeCommand = '''
FINAL_BEDPE_FILE="${{FINAL_NMSRT_BAM_PREFIX}}.bedpe.gz"
FINAL_TA_FILE="${{FINAL_BAM_PREFIX}}.final.tagAlign.gz"

bedtools bamtobed -bedpe -mate1 -i ${{FINAL_NMSRT_BAM_FILE}} | \\
  gzip -nc > ${{FINAL_BEDPE_FILE}}
zcat ${{FINAL_BEDPE_FILE}} | \\
  awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",$1,$2,$3,$9,$4,$5,$6,$10}}' | \\
  gzip -nc > ${{FINAL_TA_FILE}}

# link ${{FINAL_TA_FILE}} to ${{ipDir}}
ln -sf ${{FINAL_TA_FILE}} ${{ipDir}}

'''

step2bBothCommand = '''
#### for both PE and SE samples
# Trim R1 fastq to {trim}bp
echo "Starting trimming fastq..."

TRIM_OFPREFIX="${{prefix}}.trim"
TRIMMED_FASTQ_R1="${{TRIM_OFPREFIX}}.fastq.gz"

python $(which trimfastq.py) ${{fastqR1}} {trim} | gzip -nc > ${{TRIMMED_FASTQ_R1}}

# Align $TRIMMED_FASTQ_R1 (not paired) with bowtie2 (step 1a SE) and use it for filtering
# step (1b) and then get $FILT_BAM_FILE (not the deduped $FINAL_BAM_FILE), which is
# filtered but not deduped.
echo "Start trimmed-reads mapping..."

TRIMLOG="${{TRIM_OFPREFIX}}.align.log"
TRIM_RAW_BAM_FILE="${{TRIM_OFPREFIX}}.bam"
TRIM_FILT_BAM_PREFIX="${{TRIM_OFPREFIX}}.filt.srt"
TRIM_FILT_BAM_FILE="${{TRIM_FILT_BAM_PREFIX}}.bam"
MAPQ_THRESH={quality}

zcat -f ${{TRIMMED_FASTQ_R1}} | bowtie2 --mm -x ${{bwt2_idx}} --threads {threadNum} -U - 2> ${{TRIMLOG}} | \\
  samtools view -Su /dev/stdin | \\
  samtools sort -T ${{TRIM_OFPREFIX}} -O bam -o ${{TRIM_RAW_BAM_FILE}}
samtools view -F 1804 -q ${{MAPQ_THRESH}} -b ${{TRIM_RAW_BAM_FILE}} -o ${{TRIM_FILT_BAM_FILE}}

# =================================
# make tagAlign for filtered (but not deduped) BAM
# and subsample it for cross-correlation analysis
# ================================
echo "Make tagAlign for trim-filtered (but not deduped) BAM..."

TA_FILE="${{TRIM_FILT_BAM_PREFIX}}.tagAlign.gz"
NREADS=15000000
SUBSAMPLED_TA_FILE="${{TRIM_OFPREFIX}}.filt.sample.$((NREADS / 1000000)).tagAlign.gz"

bedtools bamtobed -i ${{TRIM_FILT_BAM_FILE}} | \\
  awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}' | \\
  gzip -nc > ${{TA_FILE}}
zcat -f ${{TA_FILE}} | grep -v “chrM” | \
  shuf -n ${{NREADS}} --random-source=<(openssl enc -aes-256-ctr -pass \\
  pass:$(zcat -f ${{TA_FILE}} | wc -c) -nosalt </dev/zero 2>/dev/null) | \\
  gzip -nc > ${{SUBSAMPLED_TA_FILE}}

### cross-correlation analysis
# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak
# <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab>
# relPhantomPeakCoef <tab> QualityTag

### temporary block
CC_SCORES_FILE="${{SUBSAMPLED_TA_FILE}}.cc.qc"
CC_PLOT_FILE="${{SUBSAMPLED_TA_FILE}}.cc.plot.pdf"
NTHREADS={threadNum}

Rscript $(which run_spp.R) -c=${{SUBSAMPLED_TA_FILE}} -p=${{NTHREADS}} -filtchr=chrM \\
  -savp=${{CC_PLOT_FILE}} -out=${{CC_SCORES_FILE}}
sed -r 's/,[^\\t]+//g' ${{CC_SCORES_FILE}} > temp
mv temp ${{CC_SCORES_FILE}}

'''

step2cSeCommand = '''
PR_PREFIX="${{OFPREFIX}}.filt.nodup"
PR1_TA_FILE="${{PR_PREFIX}}.pr1.tagAlign.gz"
PR2_TA_FILE="${{PR_PREFIX}}.pr2.tagAlign.gz"

# Get total number of read pairs
nlines=$( zcat -f ${{FINAL_TA_FILE}} | wc -l )
nlines=$(( (nlines + 1) / 2 ))

# Shuffle and split BED file into 2 equal parts

zcat -f ${{FINAL_TA_FILE}} | shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f \\
  ${{FINAL_TA_FILE}} | wc -c) -nosalt </dev/zero 2>/dev/null) | split -d -l ${{nlines}} - ${{PR_PREFIX}}

# Will produce ${{PR_PREFIX}}00 and ${{PR_PREFIX}}01
# Convert reads into standard tagAlign file
gzip -nc "${{PR_PREFIX}}00" > ${{PR1_TA_FILE}}
rm "${{PR_PREFIX}}00"
gzip -nc "${{PR_PREFIX}}01" > ${{PR2_TA_FILE}}
rm "${{PR_PREFIX}}01"

# link ${{PR1_TA_FILE}} and ${{PR2_TA_FILE}} to ${{ipDir}}
ln -sf ${{PR1_TA_FILE}} ${{ipDir}}
ln -sf ${{PR2_TA_FILE}} ${{ipDir}}

'''

step2cPeCommand = '''
PR_PREFIX="${{OFPREFIX}}.filt.nodup"
PR1_TA_FILE="${{PR_PREFIX}}.pr1.tagAlign.gz"
PR2_TA_FILE="${{PR_PREFIX}}.pr2.tagAlign.gz"
joined="${{FINAL_BAM_PREFIX}}.final.tagAlign.temp.bedpe"

# Make temporary fake BEDPE file from FINAL_TA_FILE
zcat -f ${{FINAL_TA_FILE}} | sed 'N;s/\\n/\\t/' > $joined

# Get total number of read pairs
nlines=$( zcat -f ${{joined}} | wc -l )
nlines=$(( (nlines + 1) / 2 ))

# Shuffle and split BEDPE file into 2 equal parts
zcat -f ${{joined}} | shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f \\
  ${{FINAL_TA_FILE}} | wc -c) -nosalt </dev/zero 2>/dev/null) | split -d -l ${{nlines}} - ${{PR_PREFIX}}

# Will produce ${{PR_PREFIX}}00 and ${{PR_PREFIX}}01
# Convert fake BEDPE into standard tagAlign file
awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",
  $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' "${{PR_PREFIX}}00" | gzip -nc > ${{PR1_TA_FILE}}
rm "${{PR_PREFIX}}00"
awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",
  $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' "${{PR_PREFIX}}01" | gzip -nc > ${{PR2_TA_FILE}}
rm "${{PR_PREFIX}}01"
rm -f ${{joined}}

# link ${{PR1_TA_FILE}} and ${{PR2_TA_FILE}} to ${{ipDir}}
ln -sf ${{PR1_TA_FILE}} ${{ipDir}}
ln -sf ${{PR2_TA_FILE}} ${{ipDir}}

'''

step2dBothCommand = '''
POOLED_TA_FILE="${{DATASET_PREFIX}}.pooled.tagAlign.gz"
FINAL_TA_FILES=""
for i in `find ./ -maxdepth 1 -type l -name "*.final.tagAlign.gz"|grep -v 'pooled'|sort`;
do
  FINAL_TA_FILES=`echo "${{FINAL_TA_FILES}}"" ""${{i}}"`
done

zcat -f ${{FINAL_TA_FILES}} | gzip -nc > ${{POOLED_TA_FILE}}

# ========================
# Create pooled pseudoreplicates
# =======================
# PR_PREFIX="${{OFPREFIX}}.filt.nodup"
# PR1_TA_FILE="${{PR_PREFIX}}.pr1.tagAlign.gz"
# PR2_TA_FILE="${{PR_PREFIX}}.pr2.tagAlign.gz"

PPR1_TA_FILE="${{DATASET_PREFIX}}.pooled.pr1.tagAlign.gz"
PPR2_TA_FILE="${{DATASET_PREFIX}}.pooled.pr2.tagAlign.gz"

PPR1_TA_FILES=""
for i in `find ./ -maxdepth 1 -type l -name "*.pr1.tagAlign.gz"|grep -v 'pooled'|sort`;
do
  PPR1_TA_FILES=`echo "${{PPR1_TA_FILES}}"" ""${{i}}"`
done

PPR2_TA_FILES=""
for i in `find ./ -maxdepth 1 -type l -name "*.pr2.tagAlign.gz"|grep -v 'pooled'|sort`;
do
  PPR2_TA_FILES=`echo "${{PPR2_TA_FILES}}"" ""${{i}}"`
done

zcat -f ${{PPR1_TA_FILES}} | gzip -nc > ${{PPR1_TA_FILE}}
zcat -f ${{PPR2_TA_FILES}} | gzip -nc > ${{PPR2_TA_FILE}}

'''

step2fBothCommand = '''
NTH={threadNum}
BLACKLIST={blacklist}
MAPQ_THRESH={quality}
JSD_LOG="${{DATASET_PREFIX}}.JSD.log"
JSD_PLOT="${{DATASET_PREFIX}}.JSD.png"

# BAMs are blacklist-filtered first for each replicate and control
# FINAL_BAM_PREFIX="${{DATASET_PREFIX}}.filt.srt.nodup"
# FINAL_BAM_FILE="${{FINAL_BAM_PREFIX}}.bam" # To be stored

echo "Filtering bam with blacklist bed file..."
NODUP_BFILT_BAMS=""
REPS=""
COUNTER=1
for i in `find ./ -maxdepth 1 -type l -name "*.filt.srt.nodup.bam"|grep -v 'pooled'|sort`;
do
  BTILT_PREFIX=${{i%%.bam}}
  NODUP_BFILT_BAM=${{BTILT_PREFIX}}.bfilt.bam
  FINAL_BAM_INDEX_FILE=${{BTILT_PREFIX}}.bfilt.bai
  bedtools intersect -nonamecheck -v -abam ${{i}} -b ${{BLACKLIST}} > ${{NODUP_BFILT_BAM}}
  samtools index ${{NODUP_BFILT_BAM}} ${{FINAL_BAM_INDEX_FILE}}
  NODUP_BFILT_BAMS=`echo "${{NODUP_BFILT_BAMS}}"" ""${{NODUP_BFILT_BAM}}"`
  REPS=`echo "${{REPS}}"" ""rep${{COUNTER}}"`
  COUNTER=$[$COUNTER +1]
done

echo "Plotting Jensen-Shannon distance"

C_ALL=en_US.UTF-8 LANG=en_US.UTF-8 plotFingerprint -b ${{NODUP_BFILT_BAMS}} \\
  --labels ${{REPS}} --outQualityMetrics ${{JSD_LOG}} \\
  --minMappingQuality ${{MAPQ_THRESH}} -T "Fingerprints of different samples" \\
  --numberOfProcessors ${{NTH}} --plotFile ${{JSD_PLOT}}

'''

step2gBothCommand = '''
# we don't use plot directly generated from picard
# we process picard’s text output and make a plot
# FINAL_BAM_PREFIX="${{OFPREFIX}}.filt.srt.nodup"
# FINAL_BAM_FILE="${{FINAL_BAM_PREFIX}}.bam" # To be stored

GCBIAS="/opt/picard/2.21.1/picard.jar CollectGcBiasMetrics"
GC_BIAS_LOG="${{FINAL_BAM_PREFIX}}.gc_bias.txt"
GC_BIAS_PLOT="${{FINAL_BAM_PREFIX}}.gc_bias.pdf"
SUMMARY="${{FINAL_BAM_PREFIX}}.gc_bias.summary.txt"
fasta={fasta}

java -Xmx6G -XX:ParallelGCThreads=1 -jar \\
  ${{GCBIAS}} R=${{fasta}} I=${{FINAL_BAM_FILE}} O=${{GC_BIAS_LOG}} \\
  USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE \\
  VERBOSITY=ERROR QUIET=TRUE \\
  ASSUME_SORTED=FALSE \\
  CHART=${{GC_BIAS_PLOT}} S=${{SUMMARY}}
# use ${{GC_BIAS_LOG}} into the following pyhton script
# data_file: ${{GC_BIAS_LOG}}
# prefix: any good prefix for output file name

source activate py3
plotGC.py -data ${{GC_BIAS_LOG}} -prefix ${{FINAL_BAM_PREFIX}}

'''

# ================================
#sbatch main template: 
#step 1 to 2c, 2g
# ================================
sbatchS1To2cgTemplate = '''#!/bin/bash
#SBATCH --job-name=ChipAlign_{baseExpName}    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@cho.org          # Mail user
#SBATCH -n {threadNum}                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem={memory}                      # Amount of memory in GB
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output={s1T2cLogFilePath}   # Standard output and error log

### this pipeline mainly based on ENCODE3 pipeline:
###https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c

# ================================
# load required module
# ================================
echo -e "load modules..."
module load picard/2.21.1
module load bedtools/2.29.0
module load Python/3.6.5
module load Anaconda2/4.2.0

# ===========step 0a==============
# public bariables
# ================================
echo -e "\\nPreparing basic variables..."

BASE="{basepath}"
mainAlignDir="{mainAlignDir}"
ipDir="${{mainAlignDir}}/{exp}/{ip}"
prefix="${{ipDir}}/{rep}/{baseExpName}"
fastqR1="{fastqR1}"
fastqR2="{fastqR2}"
bwt2_idx="{genomeIndex}"
log="${{prefix}}.align.log"
flagstat_qc="${{prefix}}.flagstat.qc"

# ===========step 0b==============
# cat multiple run fastqs if necessary 
# ================================
{catCommand}

# ===========step 1a==============
# Read alignment {pairedStatus} (bowtie2 aligner)
# ================================
# output:
# ● BAM file $bam
# ● mapping stats from flagstat (SAMstats) $flagstat_qc
echo "Step 1a starting..."
echo -e "Start read alignment (bowtie2 aligner)..."

{step1aCommand}
echo "Step 1a done!"

# =============step 1b==============
# Post-alignment filtering #
# Remove unmapped, mate unmapped
# not primary alignment, reads failing platform
# Remove low MAPQ reads
# Only keep properly paired reads
# Obtain name sorted BAM file
# ================================
# output:
# ● Filtered deduped position sorted BAM and index file ${{FINAL_BAM_FILE}}
# ${{FINAL_BAM_INDEX_FILE}}
# ● Filtered deduped name sorted BAM file ${{FINAL_NMSRT_BAM_FILE}}
# ● Flagstat Metric for filtered BAM file ${{FINAL_BAM_FILE_MAPSTATS}}
# ● Duplication metrics from MarkDuplicates ${{DUP_FILE_QC}}
# ● Library complexity measures ${{PBC_FILE_QC}}
echo "Step 1b starting..."
echo "Start post-alignment filtering..."

{step1bCommand}
echo "Step 1b done!"

# =======step 2a=================
# Convert SE & PE BAM to tagAlign (BED 3+3 format) #
# Create tagAlign file
# Create BEDPE file (for PE dataset)
# =================================
# output:
# ● tagAlign file ${{FINAL_TA_FILE}}
echo "Create tagAlign (& BEDPE) file..."
echo "Step 2a starting..."

{step2aCommand}
echo "Step 2a done!"

# ===========step 2b==============
#Calculate Cross-correlation QC scores
# =======================================
# output:
# ● outFile containing NSC/RSC results in tab-delimited file of 11 columns (same file can
# be appended to from multiple runs) ${{CC_SCORES_FILE}}
# ● cross-correlation plot ${{CC_PLOT_FILE}}
echo "Step 2b starting..."
echo "Calculating Cross-correlation QC scores..."

{step2bCommand}
echo "Step 2b done!"

# ==========step 2c==================
# Create pseudoReplicates
# input TagAlign file $FINAL_TA_FILE
# =====================================
# output:
# 2 pseudoreplicate virtual SE tagAlign files ${{PR1_TA_FILE}} ${{PR2_TA_FILE}}
echo "Step 2c starting..."
echo "Creating pseudoReplicates..."

{step2cCommand}
echo "Step 2c done!"

echo "Skip step to 2g..."
# ========step 2g=============
# Calculate GC bias
# =======================
# output:
# ● GC bias Plot ${{GC_BIAS_PLOT}}
# ● GC bias log ${{GC_BIAS_LOG}}
echo "Step 2g starting..."
{step2gCommand}
echo "Step 2g done!"

'''

# ================================
#sbatch main template: 
#step 2d to 3
# ================================

sbatchS2dTo2fTemplate = '''#!/bin/bash
#SBATCH --job-name=ChipPoolPeak_{baseExpName}    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@cho.org          # Mail user
#SBATCH -n 1                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem={memory}                      # Amount of memory in GB
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output={s2dTo2fLogFilePath}   # Standard output and error log

### this pipeline mainly based on ENCODE3 pipeline:
###https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c

# ================================
# load required module
# ================================
echo -e "load modules..."
module load picard/2.21.1
module load bedtools/2.29.0
module load Python/3.6.5
module load Anaconda2/4.2.0

# ===========step 0==============
# public bariables
# ================================
echo -e "\\nPreparing basic variables..."

BASE="{basepath}"
mainAlignDir="{mainAlignDir}"
ipDir="${{mainAlignDir}}/{exp}/{ip}"
DATASET_PREFIX="${{ipDir}}/{baseExpName}"
OFPREFIX="${{ipDir}}/{rep}/{baseExpName}"

echo "Change wd directory to ${{ipDir}}"
cd ${{ipDir}}

# =========step 2d===========
# Create pooled datasets
# =======================
# output:
# ● Pooled tagAlign file ${{POOLED_TA_FILE}}
# ● 2 pooled-pseudoreplicate tagAlign files
echo "Step 2d starting..."
{step2dCommand}
echo "Step 2d done!"

# ========step 2f============
# Calculate Jensen-Shannon distance (JSD)
# =======================
# output:
# ● JSD Plot ${{JSD_PLOT}}
# ● JSD log ${{JSD_LOG}}
echo "Step 2f starting..."
{step2fCommand}
echo "Step 2f done!"

'''

# generate sbatch scripts and runBash script
runSbatchScript = os.path.join(bashDir, "runEncodeChipBowtie2Qc.sbatch")
with open(runSbatchScript, 'w') as sbatchO:
    sbatchO.write('#!/bin/sh\n')
    sbatchO.write('BASE={0}\n\n'.format(bashDir))
    sbatchO.write('#sbatch for step 1 to step 2c and step 2g\n')
    # expFqDict: exp->IP->rep->[fq1,fq2]
    for exp in sorted(expFqDict.keys()):
        for ip in sorted(expFqDict[exp].keys()):
            ipDict = expFqDict[exp][ip]
            for rep in sorted(ipDict.keys()):
                pairedNumList = sorted(ipDict[rep].keys())
                baseExpName = '_'.join([exp, ip, rep])
                baseAlignDir = os.path.join(mainAlignDir, exp, ip, rep)
                os.makedirs(baseAlignDir, exist_ok=True)
                ##cat fastq files for multiple runs
                pairedNum1 = pairedNumList[0]
                runFileList = sorted(ipDict[rep][pairedNum1])
                catCommand = ''
                if len(pairedNumList) > 1:
                    if len(runFileList) > 1:
                        catFastqR1 = '${prefix}' + '_1.{fqExt}'.format(**vars())
                        catFastqR2 = '${prefix}' + '_2.{fqExt}'.format(**vars())
                        fastqR1Runs = ' '.join(list(map(lambda x: '$BASE/'+x, sorted(ipDict[rep]['1']))))
                        fastqR2Runs = ' '.join(list(map(lambda x: '$BASE/'+x, sorted(ipDict[rep]['2']))))
                        catCommand = 'echo "cat fastqs..."\ncat {fastqR1Runs} | \\ >{catFastqR1}\n\
                            cat {fastqR2Runs} | \\ > {catFastqR2}\n'.format(**vars())
                        ipDict[rep]['1'] = catFastqR1
                        ipDict[rep]['2'] = catFastqR2
                    else:
                        ipDict[rep]['1'] = '$BASE/' + ipDict[rep]['1'][0]
                        ipDict[rep]['2'] = '$BASE/' + ipDict[rep]['2'][0]
                else:
                    if len(runFileList) > 1:
                        catFastq =  '${prefix}' + '.fastq'
                        fastqRuns = ' '.join(list(map(lambda x: '$BASE/'+x, sorted(ipDict[rep][pairedNum1]))))
                        catCommand = 'echo "cat fastqs..."\ncat {fastqRuns} > {catFastq}\n'.format(**vars())
                        ipDict[rep][pairedNum1] = catFastq
                    else:
                        ipDict[rep][pairedNum1] = '$BASE/' + ipDict[rep][pairedNum1][0]
                ## generate sbatch script for step1 to step2c and step 2g
                sbatchS1To2cgScript = os.path.join(bashDir, baseExpName + '.s1To2cg.sh')
                with open (sbatchS1To2cgScript, 'w') as sbatchS1To2cg:
                    step2bCommand = step2bBothCommand.format(**vars())
                    step2gCommand = step2gBothCommand.format(**vars())
                    if len(pairedNumList) > 1:
                        fastqR1 = ipDict[rep]['1']
                        fastqR2 = ipDict[rep]['2']
                        pairedStatus = "Paired-end"
                        step1aCommand = setp1aPeCommand.format(**vars())
                        step1bCommand = step1bPeCommand.format(**vars())
                        step2aCommand = step2aPeCommand.format(**vars())
                        step2cCommand = step2cPeCommand.format(**vars())
                    else:
                        fastqR1 = ipDict[rep][pairedNum1]
                        fastqR2 = ""
                        pairedStatus = "Single-end"
                        step1aCommand = step1aSeCommand.format(**vars())
                        step1bCommand = step1bSeCommand.format(**vars())
                        step2aCommand = step2aSeCommand.format(**vars())
                        step2cCommand = step2cSeCommand.format(**vars())
                    s1T2cLogFileName = '_'.join(['runEncodeChipBowtie2Qc', baseExpName, '.s1To2cg.log'])
                    s1T2cLogFilePath = os.path.join(basepath, 'runBash', 'log', s1T2cLogFileName)
                    sbatchS1To2cgCont = sbatchS1To2cgTemplate.format(**vars())
                    sbatchS1To2cg.write(sbatchS1To2cgCont)
                    os.chmod(sbatchS1To2cgScript, 0o755)
                    sbatchO.write('sbatch $BASE/{baseExpName}.s1To2cg.sh\n'.format(**vars()))
    ## generate sbatch script for step1 to step2d and step 2f
    sbatchO.write('\n#sbatch for step 2d to step 2f\n')
    for exp in sorted(expFqDict.keys()):
        for ip in sorted(expFqDict[exp].keys()):
            ipDict = expFqDict[exp][ip]
            baseExpName = '_'.join([exp, ip])
            sbatchS2dTo2fScript = os.path.join(bashDir, baseExpName + '.s2dTo2f.sh')
            with open (sbatchS2dTo2fScript, 'w') as sbatchS2dTo2f:
                step2dCommand = step2dBothCommand.format(**vars())
                step2fCommand = step2fBothCommand.format(**vars())
                s2dTo2fLogFileName = '_'.join(['runEncodeChipBowtie2Qc', baseExpName, '.s2dTo2f.log'])
                s2dTo2fLogFilePath = os.path.join(basepath, 'runBash', 'log', s2dTo2fLogFileName)
                sbatchS2dTo2fCont = sbatchS2dTo2fTemplate.format(**vars())
                sbatchS2dTo2f.write(sbatchS2dTo2fCont)
                os.chmod(sbatchS1To2cgScript, 0o755)
                sbatchO.write('sbatch $BASE/{baseExpName}.s2dTo2f.sh\n'.format(**vars()))
os.chmod(runSbatchScript, 0o755)
