### Here is another very good source for Rainfall & other genomic plots in R
## @ https://www.bioconductor.org/packages/3.3/bioc/vignettes/gtrellis/inst/doc/gtrellis.html
library(circlize)
par(mar = c(1, 1, 1, 1))
load(paste(system.file(package = "circlize"), "/extdata/DMR.RData", sep=""))
pdf("CicularDiffMetPlot.pdf", width = 10, height = 10)
circos.initializeWithIdeogram(plotType = c("ideogram","labels"))
bed_list = list(DMR_hyper, DMR_hypo)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))
#circos.genomicDensity(bed_list[[1]], col = c("#FF000080"), track.height = 0.1)
#circos.genomicDensity(bed_list[[2]], col = c("#0000FF80"), track.height = 0.1)
dev.off()
