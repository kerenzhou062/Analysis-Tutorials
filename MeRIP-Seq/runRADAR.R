library(RADAR)

bam_dir <- "\\\\isi-dcnl\\user_data\\JJChen_Grp\\zhoukr\\project\\brandontan\\fto\\m6a_seq\\inhouse_data\\bams\\shaFTO"
samplenames <- c('shaFTO1', 'shaFTO2', 'shNS1', 'shNS2')
samplenames
radar <- countReads(samplenames = samplenames,
                    gtf = "\\\\isi-dcnl\\user_data\\JJChen_Grp\\zhoukr\\public\\genome\\annotation\\hg38\\v33\\gencode.v33.annotation.gtf",
                    bamFolder = bam_dir,
                    modification = "m6A",
                    fragmentLength = 150,
                    outputDir = "test",
                    threads = 8,
                    binSize = 50
)

radar_bak <- radar


summary(radar)
radar <- normalizeLibrary(radar)
sizeFactors(radar)
#radar <- adjustExprLevel(radar)

radar <- adjustExprLevel(radar)
variable(radar) <- data.frame( Group = c( "shaFTO", "shaFTO", "shNS", "shNS" ) )
radar <- filterBins(radar,minCountsCutOff = 5)

#radar <- diffIP(radar)
radar <- diffIP_parallel(radar, thread = 4)

#top_bins <- extractIP(radar,filtered = T)[order(rowMeans(extractIP(radar,filtered = T) ),decreasing = T)[1:1000],]

#plotPCAfromMatrix(top_bins,group = unlist(variable(radar)))

radar <- reportResult(radar, cutoff = 0.1, Beta_cutoff = 0.5,threads = 8)
result <- results(radar)
resultFile <- file.path('.\\', paste("radar", ".sig.txt", sep=""))
output.file <- file(resultFile, "wb")
write.table(as.data.frame(result), sep="\t", eol = "\n", 
            quote = FALSE, row.names=FALSE, file=output.file)
close(output.file)

#plotHeatMap(radar)

