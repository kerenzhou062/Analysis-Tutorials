# A Short Tutorial for longRNA-seq Data (edit from ENCODE pipeline)

[Ke-Ren Chow](https://github.com/KR-Chow/)

* * *

## Table of Contents
This pipeline including quality cheking, reads mapping and estimation of gene expressions.

* [Introduction](#intro)
* [Quality check](#quality)
* [Mapping](#mapping)
* [Estimation of Gene Expression](#expression)
* [Differential Gene Expression](#peakCalling)
* [Reference](#ref)
* [Licence](#licence)
* [Acknowledgments](#acknowledgments)
* * *

## <a name="intro"></a> Introduction
RNA-Seq is used to analyze the continuously changing cellular transcriptome. 
Specifically, RNA-Seq facilitates the ability to look at alternative gene spliced transcripts, 
post-transcriptional modifications, gene fusion, mutations/SNPs and changes in gene expression over time, 
or differences in gene expression in different groups or treatments. In addition to mRNA transcripts, 
RNA-Seq can look at different populations of RNA to include total RNA, small RNA, such as miRNA, tRNA, 
and ribosomal profiling. RNA-Seq can also be used to determine exon/intron boundaries and verify or amend previously
annotated 5' and 3' gene boundaries. Recent advances in RNA-seq include single cell sequencing and in situ sequencing of fixed tissue.

## <a name="quality"></a> Quality Check

Before mapping reads to the genome, we should check the qualities of sequencing results. 

The qualities of sequencing reads were checked using [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

```
fastqc fastq1.fastq fastq2.fastq -o outputDirectory
```

## <a name="mapping"></a> Reads Mapping 

After confirming reads passed quality check, we align sequencing reads to the genome using [STAR](https://github.com/alexdobin/STAR) (ENCODE recommended) or [HISAT2](http://bowtie-bio.sourceforge.net/index.shtml).

### Mapping using STAR
* installation
```bash
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.1a.tar.gz
tar -xzf 2.7.1a.tar.gz
cd STAR-2.7.1a

# Compile
cd STAR/source
make STAR
# for easy use, add bin/ to your PATH
```

* building index of reference genome
```bash
# downloading dna index fasta file
nohup wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz &

# download gft annotation file
nohup wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz &

mkdir STAR_index && cd STAR_index
STAR --runMode genomeGenerate --genomeDir STAR_index/ --genomeFastaFiles hg38.fa --sjdbGTFfile gencode.v31.annotation.gtf --sjdbOverhang 199

# --sjdbOverhang = reads length-1
```

* basic usage of STAR
```bash
STAR
--runThreadN NumberOfThreads
--genomeDir /path/to/genomeDir
--readFilesIn /path/to/read1 [/path/to/read2]
--outSAMtype BAM SortedByCoordinate

# example
STAR --runThreadN 20 --genomeDir STAR_index/ --readFilesCommand zcat \
--readFilesIn reads_1.fq.gz reads_2.fq.gz --outSAMtype BAM SortedByCoordinate

# unsorted or sorted bam file
--outSAMtype BAM Unsorted, or
--outSAMtype BAM SortedByCoordinate, or
--outSAMtype BAM Unsorted SortedByCoordinate, or

```
* ENCODE pipline: [STAR_RSEM.sh](https://github.com/KR-Chow/Analysis-Tutorials/blob/master/scripts/STAR_RSEM.sh)

### Mapping using HISAT2
* installation
```bash
# Get latest STAR source from releases
wget http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-Linux_x86_64.zip
tar -xzf hisat2-2.1.0-Linux_x86_64.zip
cd hisat2-2.1.0

# for easy use, add bin/ to your PATH
```

* building index of reference genome
```bash
# Extract splice-site and exon information from gene annotation files
extract_splice_sites.py gencode.v31.annotation.gtf > gencodeV31.ss
extract_exons.py gencode.v31.annotation.gtf > gencodeV31.exon

# Building a HISAT2 index
hisat2-build --ss gencodeV31.ss --exon gencodeV31.exon hg38.fa hg38Transcrits
```

* align reads to reference genome
```bash
# for paired-end reads
hisat2 -p 3 --rna-strandness R --dta -x hg38Transcrits -q -1 reads_1.fq -2 reads_2.fq -S readsMapped.sam

# for single reads
hisat2 -p 3 --rna-strandness R --dta -x hg38Transcrits -q reads.fq -S readsMapped.sam

# sam to index bam by using samtools
sort --threads 8 -O bam -o reads.sorted.bam reads.sam
samtools index -b reads.sorted.bam

# stringtie: assemble transcriptome
stringtie -p 8 --rf -G gencode.v31.annotation.gtf -o readsMapped.gtf readsMapped.sorted.bam
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
