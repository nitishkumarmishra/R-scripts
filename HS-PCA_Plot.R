library(readxl)
library(ggplot2)
library(dplyr)

setwd("C:/Users/nitis/Desktop/Radha HS")
Data <- as.data.frame(read_excel("12-HS group FDR and AUC removed dupliates Feb 18 2021 to Nitish Heatmap.xlsx", sheet = 1))
rownames(Data) <- Data$TargetID...2
meth <- Data[grepl("AVG_Beta", colnames(Data))]
Group1 <- grepl("HS", colnames(meth))
Group1 <- meth[,Group1]

Group2 <- grepl("CTR", colnames(meth))
Group2 <- meth[,Group2]
sample_set <- cbind(Group1, Group2)
## PCA Plot
pca <- prcomp(t(sample_set), scale. = TRUE)
pd <- ifelse(grepl("HS", colnames(sample_set)),"HS", "Control")
color <- as.factor(pd)
PCi<-data.frame(pca$x,Sample=color)
## 2D plot
ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("HS"="red4","Control"="blue4"))+ #your colors here
  theme_classic()
ggsave("2D_PCA_selected.pdf", dpi = 600, width = 12, height = 8)
################### 3D PCA plot ######################
groups <- levels(color)
pca$pcolor <- pd
pca$pcolor <- ifelse(grepl("HS", pca$pcolor), "red4", "blue4")
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 60)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("topright", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
dev.print(pdf, '3D_PCA_All.pdf', width = 10, height = 10)
########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 60)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("topright", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100),labels = gsub(".AVG_Beta", "",rownames(pca$x)), cex= 0.6, col = pca$pcolor)
dev.print(pdf, '3D_PCA_withXY_withName.pdf', width = 12, height = 12)

############## Volcano plot ##########
Data$Methylation.diff <- (Data$`% Methylation  Cases` - Data$`% Methylation Control`)/100
x.cut=0.1;y.cut=0.01
Significance <- ifelse(Data$Methylation.diff >= x.cut & Data$`FDR p-Val` < y.cut, "Hypermethylated", ifelse(Data$Methylation.diff <= -x.cut & Data$`FDR p-Val` < y.cut, "Hypomethylated", "Not significant"))
Gene_Type <- "methylation"
ggplot(Data, aes(x = Methylation.diff, y = -log10(`FDR p-Val`))) +
  geom_point(aes(color = Significance), size = 1.2) +
  scale_color_manual(values = c( "blue4","red4","gray41" )) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  #geom_text_repel(data = subset(meth, ( FDR.p.Val< 0.00001 & abs(Methylation.diff) >=12)),
  #               aes(label = TargetID)) +
  
  theme(legend.position="right")+
  ggtitle("Volcano plot")+
  xlab("Delta Beta Value in Percentage") + ylab("-log10(padj)") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 

ggsave(filename = paste0(Gene_Type,"_VolcanoPlot.pdf"), width = 12, height = 8, dpi = 800)
################## Heatmap plot ############

library(tidyheatmap); library(dplyr); library(tidyverse); library(matrixStats)

Data$Status <- ifelse(Data$Methylation.diff >0 , "Hyper", "Hypo")
#colnames(Data) <- gsub(".AVG_Beta", "", colnames(Data))
Data$adj.P.Val <- Data$`FDR p-Val`
Data$deltaBeta <- Data$Methylation.diff
Data$feature <- sapply(strsplit(Data$UCSC_REFGENE_GROUP...16,";", 1), '[[', 1)

group <- pd
####################################
Num_Gene <- 200
Exp <- Data %>%
  slice_max(abs(Methylation.diff), n=Num_Gene) %>%
  mutate(log2BH= -log2(adj.P.Val)) %>%
  mutate(Feature = feature, "Delta Beta"=deltaBeta, "log2(BH)"=log2BH)%>%
  pivot_longer(grep("AVG_Beta", names(Data)),names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(group, Num_Gene))#rep(condition, 7); its seven gene so seventies of condition

ann_colors <- list(Group = c(HS = "red", Control = "blue"),
                   Status=c(Hyper="red", Hypo="blue"),
                   #CGI=c(island="blue", opensea="red", shelf="green", shore="darkcyan"),
                   #logFC= c("blue","gray","red"),
                   "log2(BH)"=c("blue","gray","red"),
                   "Delta Beta"=c("blue","gray","red"),
                   Feature=c(Body="darkgreen",TSS1500="red","5'UTR"="orange","3'UTR"= "magenta4",TSS200="blue","1stExon"="blueviolet"))

tidy_heatmap(Exp,
             rows = "TargetID...2",
             columns = Sample,
             values = Expression,
             scale = "none",
             annotation_col = c(Group),
             annotation_row = c(Status, "Delta Beta", "log2(BH)",Feature),
             cluster_cols = TRUE,
             clustering_method = "ward.D2",
             #clustering_distance_cols = "euclidean",
             clustering_distance_cols = "manhattan",
             cluster_rows = TRUE,
             color_legend_n = 7,
             colors = c("red","gray","green"),
             annotation_colors = ann_colors,
             fontsize = 10,
             fontsize_row = 6,
             fontsize_col = 8,
             angle_col = 315,
             height = 14,
             width = 12,
             show_rownames = FALSE, ## To show the rownames
             filename = "Heatmap_200probes.pdf"
)

###############################################
save.image("HS_PCA_Plot.RData")


