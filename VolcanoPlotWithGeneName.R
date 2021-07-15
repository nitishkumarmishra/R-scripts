library(ggplot2)
library(ggrepel)
library(TCGAbiolinks)
genes <- read.table("results.txt", header=TRUE)
rownames(genes) <- genes$Gene
x <- genes$log2FoldChange; y <- genes$padj
x.cut = 1; y.cut = 0.05#### Here we can choose Foldchange and adjusted P-value threshold
.e <- environment()
#genes$Significant <- ifelse(genes$padj < 0.05, "FDR < 0.05", "Not Sig")
##### This is another plot ########
Significant <- ifelse(genes$log2FoldChange >= x.cut & genes$padj < y.cut, "Up regulated", ifelse(genes$log2FoldChange <= -x.cut & genes$padj < y.cut, "Down regulated", "Not significant"))
p <- ggplot(genes, aes(x = log2FoldChange, y = -log10(padj)),size = 5) +
  geom_point(aes(color = Significant), size = 2) +
  scale_color_manual(values = c("red", "gray30", "green4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
################################
  theme(legend.position="right")+
  ggtitle("Volcano plot Tumor Vs. Normal")+
  theme(plot.title = element_text(size = 14, face = "bold"))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 12, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) +
  geom_text_repel(
    data = subset(genes, (padj < 0.05 & abs(log2FoldChange) >=1)),
    aes(label = Gene),size = 4, 
    show.legend = FALSE, fontface = "bold", color = "black", 
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.5, "lines")
  )
########################################################
ggsave(p, filename = "VolcanoPlotWithGeneName.pdf", width = 10, height = 10, dpi = 600)
