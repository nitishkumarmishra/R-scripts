load("C:/Users/nitish.mishra/Desktop/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffExpAnalysis.RData")
genes <- as.data.frame(resSig)
genes$Gene <- sapply(strsplit(rownames(genes),"\\|"),'[[',1)
x.cut=2;y.cut=0.01
Significant <- ifelse(genes$log2FoldChange >= x.cut & genes$padj < y.cut, "Up regulated", ifelse(genes$log2FoldChange <= -x.cut & genes$padj < y.cut, "Down regulated", "Not significant"))
p <- ggplot(genes, aes(x = log2FoldChange, y = -log10(padj)),size = 4) +
  geom_point(aes(color = Significant), size = 2) +
  scale_color_manual(values = c("red", "gray20", "green4")) +
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
    ### used p value 0.00001 log2FC >=abs(6) for name in plot
    data = subset(genes, (padj < 0.00001 & abs(log2FoldChange) >=6)),
    aes(label = Gene),size = 4, 
    show.legend = FALSE, fontface = "bold", color = "black", 
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.5, "lines")
  )
ggsave(p, filename = "VolcanoPlotmRNA.pdf", width = 12, height = 8, dpi = 800)
