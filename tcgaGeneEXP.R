library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library("impute")
############ edgeR Diff Exp genes ################
edgeR.diffExp <- read.table("edgeR-FC2.0-Zero36.tsv", header = TRUE, row.names = 1, sep = "\t")
edgeR.diffExp <- edgeR.diffExp[order(rownames(edgeR.diffExp)),]
##################################################
mRNAseq <- as.matrix(mRNAseq)
mRNAseq.knn <- impute.knn(mRNAseq, k=15, rowmax = 0.20, colmax = 0.20)
mRNAseq <- mRNAseq.knn$data
#tmp <- (rowSums(mRNAseq.normal ==0))
#tmp1 <- tmp[tmp >=3]
#names(tmp1)
######Make Tumor normal sample list###############
mRNAseq <- read.table("mRNAseq_raw_counts.txt", header = TRUE, row.names = 1, sep = "\t")
sampleIDs <- gsub("\\.", "-", sampleIDs)
colnames(mRNAseq) <- sampleIDs
samplesDat <- data.frame(matrix(nrow = length(sampleIDs), ncol = 4))
for (j in 1:length(sampleIDs)) {
  tmpRow <- unlist(strsplit(sampleIDs[j], split = "-"))
  samplesDat[sampleIDs[j], ] <- tmpRow
}
sampleIDs1 <- as.character(samplesDat[, 4])
#sampleIDs1 <- substr(sampleIDs1, 1, nchar(sampleIDs1))## It don't needed as name alread 01/11 not 01A/11A
sampleIDs1 <- as.numeric(sampleIDs1)
normalSamples <- rownames(samplesDat)[sampleIDs1 < 14 & sampleIDs1 > 9]
tumorSamples <- rownames(samplesDat)[sampleIDs1 < 6]
##################################################################
mRNAseq.normal <- mRNAseq[, c(normalSamples)]
mRNAseq.tumor <- mRNAseq[, c(tumorSamples)]

##################################################################
FoldChange <- edgeR.diffExp$logFC
direction = ifelse(edgeR.diffExp$logFC > 0, "Up", "Down")
##################################################################

#ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)), rep("Normal", length(normalSamples)))), col = list(type = c("Tumor" = "red", "Normal"="blue")))
ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)))), col = list(type = c("Tumor" = "red")))
#ha2 = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)), rep("Normal", length(normalSamples))), col = list(type = c("Tumor" = "red", "Normal"="blue")), show_legend = FALSE)
column_tree = hclust(dist(t(mRNAseq.tumor.edgeR.log)), method = "ward.D2")
ht_list = Heatmap(mRNAseq.tumor.edgeR.log, name = "GeneExpression", col = colorRamp2(c(0, 6, 13), c("blue", "gray","red")),
                  cluster_columns = column_tree, show_column_names = FALSE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 4, column_title = "GeneExpression", column_title_gp = gpar(fontsize = 10), show_row_names = FALSE, row_names_max_width = unit(4, "cm"), row_title_gp = gpar(fontsize = 12))+
  Heatmap(direction, name = "Direction", col = c("Up" = "red", "Down" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(FoldChange, name = "Fold Change", col = colorRamp2(c(-5.5,-2,2,5.5), c("blue","gray","pink", "red")), column_names_gp = gpar(fontsize = 8))
pdf("GeneExpression_clust4_ward2.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Gene expression analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", legend_title_gp = gpar(fontsize = 8, fontface = "bold"), legend_label_gp = gpar(fontsize = 8))
dev.off()

# #ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)), rep("Normal", length(normalSamples)))), col = list(type = c("Tumor" = "red", "Normal"="blue")))
# ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)))), col = list(type = c("Tumor" = "red")))
# 
# #ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)))), col = list(type = c("Tumor" = "red")))
# #ha2 = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)), rep("Normal", length(normalSamples))), col = list(type = c("Tumor" = "red", "Normal"="blue")), show_legend = FALSE)
# column_tree = hclust(dist(t(mRNAseq.tumor.edgeR.log)), method = "ward.D2")
# ht_list = Heatmap(mRNAseq.tumor.edgeR.log, name = "GeneExpression", col = colorRamp2(c(0, 6, 13), c("blue", "gray","red")),
#   cluster_columns = column_tree, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 4, column_title = "GeneExpression", column_title_gp = gpar(fontsize = 10), show_row_names = TRUE, row_names_max_width = unit(4, "cm"), row_title_gp = gpar(fontsize = 12))+
#   Heatmap(direction, name = "Direction", col = c("Up" = "red", "Down" = "blue"), column_names_gp = gpar(fontsize = 8))+
#   Heatmap(FoldChange, name = "Fold Change", col = colorRamp2(c(-5.5,0,5.5), c("blue","gray","red")), column_names_gp = gpar(fontsize = 8))
# pdf("GeneExpression_clust4_ward2.pdf", width = 10, height = 10)
# draw(ht_list, newpage = FALSE, column_title = "Gene expression analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", legend_title_gp = gpar(fontsize = 8, fontface = "bold"), legend_label_gp = gpar(fontsize = 8))
# dev.off()
