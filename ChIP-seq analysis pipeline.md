# A Short Tutorial for ChIP-seq Data

[Ke-Ren Chow](https://github.com/KR-Chow/)

* * *

## Table of Contents
This pipeline including quality cheking, reads mapping and peak calling.

* [Introduction](#intro)
* [Quality check](#quality)
* [Mapping](#mapping)
* [Genomic coverage](#coverage)
* [Peak calling](#peakCalling)
* [Reference](#ref)
* [Licence](#licence)
* [Acknowledgments](#acknowledgments)
* * *

## <a name="intro"></a> Introduction
ChIP-sequencing, also known as ChIP-seq, is a method used to analyze protein interactions with DNA. 
ChIP-seq combines chromatin immunoprecipitation (ChIP) with massively parallel DNA sequencing to identify the binding sites of DNA-associated proteins. 
It can be used to map global binding sites precisely for any protein of interest. 

## <a name="quality"></a> Quality Check

Before mapping reads to the genome, we should check the qualities of sequencing results. 

The qualities of sequencing reads were checked using [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

```
fastqc fastq1.fastq fastq2.fastq -o outputDirectory
```

## <a name="mapping"></a> Reads Mapping 

After confirming reads passed quality check, we align sequencing reads to the genome using [bowtie](http://bowtie-bio.sourceforge.net/index.shtml).

```bash
#using bowtie to map sequencing
bowtie -t -p 16 -v 2 -m 1 --best --strata --sam bowtieIndex_UCSC/hg19 fastq.map.sam
samtools view -h -bS -F 4 --threads 8 fastq.map.sam fastq.map.bam

#convert sam to bam (samtools v1.3.1)
samtools sort --threads 8 -m 2G -O bam -o fastqMap.sorted.bam fastq.map.bam
samtools index -b fastq.map.sorted.bam

#delete temp files
rm -f fastq.map.sam fastq.map.bam
```

## <a name="coverage"></a> Genomic Coverage

A tiled data file (TDF) file (.tdf) is a binary file that contains genomic coverage used for faster display in IGV.

TDF can be generated using  [IGVtools](https://software.broadinstitute.org/software/igv/igvtools).

```bash
igvtools count -z 7 -w 25 -e 250 mapBam.sorted.bam mapBam.cov.tdf genome.chrom.sizes > /dev/null 2>
```

## <a name="peakCalling"></a> Peak Calling

There are 2 types of ChIP-seq peaks, narrow peaks (e.g. transcrption factor) and broad peaks (e.g. histone modification).


#### peak types of histone modifications

Peak Type | [Histone Modification](https://www.encodeproject.org/chip-seq/histone/)
----------- | ----------
Narrow      | H2AFZ, H3ac, H3K27ac, H3K4me2, H3K4me3, H3K9ac, 
Broad       | H3F3A, H3K27me3, H3K36me3, H3K4me1, H3K79me2, H3K79me3, H3K9me1, H3K9me2, H4K20me1
Exceptions  | H3K9me3

Many software are available for peak-calling.
* [MACS](http://liulab.dfci.harvard.edu/MACS/)

```bash
## ref[1]

# for replicates
samtools merge -f mergeIP.bam map1.sorted.bam map2.sorted.bam
samtools merge -f mergeinput.bam map1.sorted.bam map2.sorted.bam

# for broad peak
macs14 -t mergeIP.bam -c mergeinput.bam  -g hs -n name --nomodel \
  --shiftsize 147 --pvalue 1e-3 -B -S --call-subpeaks > macs14.log 2>&1 &

# for narrow peak
macs14 -t mergeIP.bam -c mergeinput.bam  -g hs -n name -B -S \ 
  --call-subpeaks > macs14.log 2>&1 &

```

* [MACS2](https://github.com/taoliu/MACS)

```bash
## ref[2]

# for broad peak
macs2 callpeak --broad -B --nomodel --extsize 147 -g hs \
  -t IP1.sorted.bam IP2.sorted.bam \
  -c -t input1.sorted.bam input2.sorted.bam \
  -n macs2 > shCont.macs2.log 2>&1 &

# for narrow peak
macs2 callpeak -B -g hs -t IP1.sorted.bam IP2.sorted.bam \
  -c -t input1.sorted.bam input2.sorted.bam \
  -n macs2 > shCont.macs2.log 2>&1 &

```

* [SICER](https://github.com/dariober/SICERpy)
```bash
# for histone broad peak
SICER.py -t IP.map.bam \
  -c input.map.bam -rt 0 \
  > peak.bed 2> sicer.log &
  
awk 'BEGIN{OFS="\t";FS="\t";}{if($0~/^chr/ && $8 < 0.01){print $1,$2,$3;}}' \
  peak.bed > peak.bed.q01.bed

```

* [Qeseq](https://sourceforge.net/projects/klugerlab/files/qeseq/v0.2.2/)
```bash
# for histone broad peak
bedtools bamtobed -i IP.map.bam \
  > IP.map.bed
bedtools bamtobed -i input.map.bam \
  > input.map.bed

awk '{print $1 " " (($6=="+")?$2:$3)  " " $6;}' IP.map.bed \
  > IP.map.single.3col
awk '{print $1 " " (($6=="+")?$2:$3)  " " $6;}' input.map.bed \
  > input.map.single.3col

qeseq -s 200 -v 1 -o qeseqPeak -p 0.01 -c 5 \
  IP.map.single.3col input.map.single.3col \
  > qeseq.run.log 2>&1 &

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){if($5/$6>2 && $8<1e-5){print $1,$2-1,$3}}}' \
  qeseqPeak.txt > qeseqPeak.filter.bed

```

## <a name="coverage"></a> Reference 

[1]. Identifying ChIP-seq enrichment using MACS, Nat Protoc., 2012,   
[2]. Use Model-Based Analysis of ChIP-Seq (MACS) to Analyze Short Reads Generated by Sequencing Protein–DNA Interactions in Embryonic Stem Cells.


## <a name="license"></a> License

This pipelie is developed by KR Chow - see the [LICENSE.md](LICENSE.md) file for details

## <a name="acknowledgments"></a> Acknowledgments 

* Prof. Qu
* Prof. Yang