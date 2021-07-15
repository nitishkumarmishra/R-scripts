library(ggplot2)
library(ggrepel)
library(TCGAbiolinks)
x.cut <- 1 ## log2FC cutoff
y.cut <- 0.05 ## adj P-value cutoff
filename <- c("volcano1.pdf")
show.names <- c("significant")
genes <- read.table("results.txt", header=TRUE)
rownames(genes) <- genes$Gene
x <- genes$log2FoldChange; y <- genes$padj
label = NULL
#label=c("Not Significant", "Up regulared", "Down regulated")
xlab = " Gene expression fold change (Log2)"
legend = "State"
title = "Volcano plot (Tumor Vs. Normal)"
height = 10
width = 10
highlight.color = c("orange" )
names.size = 3
dpi = 800
#genes$Significant <- ifelse(genes$padj < 0.05, "FDR < 0.05", "Not Sig")
##### This is another plot ########
genes$Significant <- ifelse(genes$log2FoldChange > 0.5 & genes$padj < 0.05, "Up regulated", ifelse(genes$log2FoldChange < 0 & genes$padj < 0.05, "Down regulated", "Not significant"))
ggplot(genes, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "black", "green")) +
  theme_bw(base_size = 16) +
  geom_text_repel(
    data = subset(genes, (padj < 0.05 & abs(log2FoldChange) >=1)),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
##########
.e <- environment()
names <- genes$Gene
threshold <- rep("1", length(genes$Gene))
sig <- y < y.cut
sig[is.na(sig)] <- FALSE
up <- x > x.cut
up[is.na(up)] <- FALSE
if (any(up & sig)) 
  threshold[up & sig] <- "2"
down <- x < (-x.cut)
down[is.na(down)] <- FALSE
if (any(down & sig)) 
  threshold[down & sig] <- "3"
df <- data.frame(x = x, y = y, threshold = threshold)

if (!is.null(names)) {
  if (show.names == "significant") {
    idx <- (up & sig) | (down & sig)
    important <- c("2", "3")
  }
}
if (is.null(label)) {
  label = c(`1` = "Not Significant", `2` = "Up regulared", 
            `3` = "Down regulated")
}
#names(label) <- as.character(1:3)
color = c("black", "red", "green4")
p <- ggplot(data = df, aes(x = x, y = -1 * log10(y), colour = threshold), 
            environment = .e) + geom_point() + ggtitle(title) + ylab(ylab) + 
  xlab(xlab) + geom_vline(aes(xintercept = -x.cut), colour = "black", 
                          linetype = "dashed") + geom_vline(aes(xintercept = x.cut), 
                                                            colour = "black", linetype = "dashed") + geom_hline(aes(yintercept = -1 * 
                                                                                                                      log10(y.cut)), colour = "black", linetype = "dashed") + 
  scale_color_manual(breaks = as.numeric(names(label)), 
                     values = color, labels = label, name = legend) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), legend.text = element_text(size = 10), 
                     axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), 
                     legend.position = "top", legend.key = element_rect(colour = "white"))
### change font size of xlab, ylab and title
# p + theme(
#   plot.title = element_text(color="red", size=14, face="bold.italic"),
#   axis.title.x = element_text(color="blue", size=14, face="bold"),
#   axis.title.y = element_text(color="#993333", size=14, face="bold")
# )
### This line print name of genes but not in red and green box

p <- p + geom_text_repel(data = subset(df, threshold %in% important), aes(label = names[idx]), size = names.size, 
                         show.legend = FALSE, fontface = "bold", color = "black", 
                         point.padding = unit(0.3, "lines"), box.padding = unit(0.5,  "lines"))


p <- p + geom_label_repel(data = subset(df, threshold %in% important), aes(label = names[idx], fill = threshold), 
                          size = names.size, show.legend = FALSE, fontface = "bold", 
                          color = "white", box.padding = unit(0.35, "lines"), 
                          point.padding = unit(0.3, "lines")) + scale_fill_manual(values = color[as.numeric(important)])

ggsave(p, filename = filename, width = width, height = height, dpi = dpi)
    
#######################I can also generate plot by using TCGABiolinks #####
library(TCGAbiolinks)
TCGAVisualize_volcano(genes$log2FoldChange,genes$padj,
                      filename = "Case3_volcanoexp.pdf",
                      x.cut = 1,
                      y.cut = 0.05,
                      names = rownames(genes),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot",
                      width = 10, height = 10, dpi = 600, show.names = "significant")
    