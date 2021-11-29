#!/usr/bin/env Rscript
suppressMessages(library('getopt'))

command =  matrix(c(
    "help",   "h", 0, "logical", "notice: sc down suggest p 1 q 0.2",
    "dbName", "d", 1, "character", "org.Hs.eg.db|org.Mm.eg.db",
    "item",   "t", 1, "character", "ENSEMBL|SYMBOL|ENTREZID",
    "pval",   "p", 1, "numeric", "pval cutoff, default is 0.05",
    "qval",   "q", 1, "numeric", "qval cutoff default is 0.2",
    "input",  "i", 1, "character", "geneList input file",
    "pdf",    "o", 1, "character", "output pdf file",
    "type",   "y", 1, "character", "BP|MF|CC",
    "text" ,  "r", 1, "character", "output the items in txt file",
    "max" ,   "m", 1 , "numeric", "max items in pdf"
), byrow=TRUE, ncol=5)

args <- getopt(command)

if(!is.null(args$help)){
    cat(paste(getopt(command, usage = T),"\n"))
    q()
}

localLib <- .libPaths()[1]
suppressMessages(library('dplyr', lib.loc=localLib)
suppressMessages(library('ChIPseeker'))
suppressMessages(library('clusterProfiler'))
suppressMessages(library(args$dbName, character.only = TRUE))

## loading packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
require(GenomicRanges)
library(ChIPpeakAnno)

dMatrix <- read.csv(file="matrix.txt", sep='\t') # exp, peak
files <- as.list(as.character(dMatrix$peak))
names(files) <- as.list(as.character(dMatrix$exp))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

##### individual
peak <- readPeakFile(files[[1]])

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

pdf('tagHeatmap.pdf', paper = "a4r")

tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
garbage <- dev.off()

pdf('avgProf.pdf', paper = "a4r")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000,
            xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency")
garbage <- dev.off()

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 1000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
pdf('peakAnno.pdf', paper = "a4r")
plotAnnoPie(peakAnno)
garbage <- dev.off()

pdf('distance2TSS.pdf', paper = "a4r")
plotDistToTSS(peakAnno,
              title="Distribution of TET1-binding loci\nrelative to TSS")
garbage <- dev.off()


##### plot together
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

pdf('peakAnno.pdf', paper = "a4r")
plotAnnoBar(peakAnnoList, title="Genomic features of TET1-binding locus")
garbage <- dev.off()

pdf('distance2TSS.pdf', paper = "a4r")
plotDistToTSS(peakAnnoList, title="Distance of TET1-binding locus relative to TSS")
garbage <- dev.off()

pdf('KEGG.pdf', width=12, height=6)
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                         fun           = "enrichKEGG",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
garbage <- dev.off()


pdf('avgProf.pdf', paper = "a4r")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=1000, facet="row")
garbage <- dev.off()

bedFile <- read.table('Kas1_control_NA_rep1_final_peaks.narrowPeak', header=FALSE)
Kas1_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3))
bedFile <- read.table('MA93_control_NA_rep1_final_peaks.narrowPeak', header=FALSE)
MA93_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3))
bedFile <- read.table('MDSL_control_NA_rep1_final_peaks.narrowPeak', header=FALSE)
MDSL_granges <- GRanges(seqnames = bedFile$V1, ranges = IRanges(start = bedFile$V2, end = bedFile$V3))

VennForBeds <- makeVennDiagram(Peaks=list(Kas1_granges, MA93_granges, MDSL_granges),
  NameOfPeaks=c("Kasumi-1", "MA93", "MDSL"),
  fill=c("orange","brown","cyan"))

plot(VennForBeds)
dev.off()



gene <- read.table(args$input,sep="\t",header=F)

# if(args$type == "up"){
#     genelist = as.character(gene[gene$V2>1,]$V1)
# }else if(args$type == "down"){
#     genelist = as.character(gene[gene$V2<1,]$V1)
# }else{
#     genelist = as.character(gene$V1)
# }

genelist <- as.character(gene$V1)
geneID <- bitr(genelist, fromType = args$item, toType = "ENTREZID", OrgDb = args$dbName)

geneEnrich <- enrichGO(
        geneID$ENTREZID,
        OrgDb = args$dbName,
        ont = args$type,
        keyType = "ENTREZID",
        pvalueCutoff = args$pval,
        pAdjustMethod = "BH",
        qvalueCutoff = args$qval,
        readable = T
    )

# exl = gsub("\\.pdf",".txt",args$out)
#h <- 5 + args$max/10
pdf(args$pdf,width=10, heigh=5)
    dotplot(geneEnrich, showCategory = args$max)
dev.off()

write.table(geneEnrich, args$text, sep = "\t")

#print("done")

