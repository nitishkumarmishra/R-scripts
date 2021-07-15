require(ggplot2)
require(ggbio)
require(biomaRt)
require(GenomeGraphs)
require(GenomicRanges)
library(BSgenome)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg19)

# get chromosome information
chr.len = seqlengths(Hsapiens)  # get chromosome lengths
chrom.length = chr.len[grep("_|M|X|Y", names(chr.len), invert = T)] # remove X,Y,M and random chromosomes if required
names(chrom.length) <- gsub("chr", "", names(chrom.length))
# create the ideogram accordingly
myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, width = chrom.length))
seqlevels(myIdeo) = names(chrom.length)
seqlengths(myIdeo) = (chrom.length)

# create object with genomic locations of methylation sites, as well as required statistics
ab_data <- read.table("Probeinfo_SAMR_0.1_1.txt", header = TRUE, sep = "\t", row.names = 1)
map_data <- ab_data
#map_data <- ab_data[c(1,5,5,4,11,12,13)] # subset data to include everything required
#colnames(map_data)[2] <- "start"
#colnames(map_data)[3] <- "end"

#map_data$CHR <- paste("chr", map_data$CHR, sep="") # ensure chromosome labels are correct

# generate the 'inner' and 'outer' objects to plot (hypo and hyper methlation)
increase <- map_data[map_data$median_diff >= 0,]
decrease <- map_data[map_data$median_diff <= 0,]
#increase <- map_data[map_data$median_diff >= 0,]
#decrease <- map_data[map_data$median_diff <= 0,]
increase$strand = "*"
decrease$strand = "*"
colnames(increase)[2] <- "chr"
colnames(decrease)[2] <- "chr"
increase <- na.omit(increase) # ensure no NA values
decrease <- na.omit(decrease) # ensure no NA values
# can set a genome wide threshold if you don't want to plot every CpG site
#decrease2 <- decrease[decrease$t.pval <= 1e-7,]
#increase2 <- increase[increase$t.pval <= 1e-7,]
increase2 <- increase
decrease2 <- decrease
# create the granges objects
g.per.inc <- with(increase2, GRanges(chr, IRanges(start, end), strand = strand))
g.per.inc = keepSeqlevels(g.per.inc, names(chrom.length))
values(g.per.inc)$id = "Hypermethylated"
values(g.per.inc)$p.val = increase2$t.pval
values(g.per.inc)$t.stat = increase2$t.stat
g.per.inc@seqinfo@seqlengths <- seqlengths(myIdeo)
#
g.per.dec <- with(decrease2, GRanges(chr, IRanges(start, end), strand = strand))
g.per.dec = keepSeqlevels(g.per.dec, names(chrom.length))
values(g.per.dec)$id = "Hypomethylated"
values(g.per.dec)$p.val = decrease2$t.pval
values(g.per.dec)$t.stat = decrease2$t.stat
g.per.dec@seqinfo@seqlengths <- seqlengths(myIdeo)

# circle plot using ggplot
pdf("Figure 1a.pdf", width = 10, height = 10)
p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", 
                              radius = 57, trackWidth = 3)

p <- p + layout_circle(c(g.per.inc, g.per.dec), geom = "point", size = 1.50,  aes(x = midpoint, y = t.stat, color = id), 
                       radius = 34, trackWidth = 50) + scale_colour_manual(values = c("red4", "blue4")) 

p$labels$colour = "" # label key accordingly
#p$labels$colour = "Differential Methylation" # label key accordingly
p + layout_circle(myIdeo, geom = "text", size = 10, face = "bold", aes(label = seqnames),
                  vjust = 0, radius = 85, trackWidth = 7) + 
  theme(title = element_text("ESR Ab methylation analysis"), 
        legend.text = element_text(colour="black", size = 14),
        legend.title = element_text(colour="black", size=14)) 
dev.off()
#END
###################################################################
