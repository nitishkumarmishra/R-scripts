library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Homo.sapiens)
keytypes(Homo.sapiens)
txs <- transcriptsBy(Homo.sapiens, 'gene', col='GENEID')
library(FDb.InfiniumMethylation.hg19)
hm450 <- get450k()
probenames <- c("cg08291996", "cg16392865", "cg00395291", "cg09310185", "cg21749424")
probes <- hm450[probenames]
result <- getNearestTSS(probes)
# head(result)
# queryHits subjectHits distance nearestGeneSymbol nearestTranscript
# cg16392865         1        3932     3415             CNPY3        uc003osy.2
# cg00395291         2        3964    10642           SLC12A7        uc003jbu.3
# cg09310185         3       35174      314             BIRC7        uc010gkc.1
# cg21749424         4       39843    17271         LINC00473        uc011egl.1
sampleIDs <- colnames(PAAD.RSEM.log2)
sampleIDs <- gsub("\\.", "-", sampleIDs)
colnames(PAAD.RSEM.log2) <- sampleIDs
samplesDat <- data.frame(matrix(nrow = length(sampleIDs), ncol = 7))
rownames(samplesDat) <- sampleIDs
for (j in 1:length(sampleIDs)) {
  tmpRow <- unlist(strsplit(sampleIDs[j], split = "-"))
  samplesDat[sampleIDs[j], ] <- tmpRow
}
sampleIDs1 <- as.character(samplesDat[, 4])
sampleIDs1 <- substr(sampleIDs1, 1, nchar(sampleIDs1) - 1) 
sampleIDs1 <- as.numeric(sampleIDs1)
normalSamples <- rownames(samplesDat)[sampleIDs1 < 20 & sampleIDs1 > 9]
tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]
#mRNAseq <- mRNAseq[, c(normalSamples, tumorSamples)]
mRNAseq <- mRNAseq[, c(tumorSamples)]
#######################################################################
############## DNA methylation heatmap ################################
ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)))), col = list(type = c("Tumor" = "red")))
ha2 = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)))), col = list(type = c("Tumor" = "red")), show_legend = FALSE)
column_tree = hclust(dist(t(mat_meth.impute.TSS.tumor)), method = "ward.D2")
ht_list = Heatmap(mat_meth.impute.TSS.tumor, name = "Methylation", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
                  cluster_columns = column_tree, show_column_names = FALSE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 4, column_title = "Methylation", column_title_gp = gpar(fontsize = 10), show_row_names = FALSE,row_title_gp = gpar(fontsize = 10))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(mat_expr.impute.TSS.tumor[, column_tree$order], name = "Expression", col = colorRamp2(c(-2.01, 0, 16.87), c("blue","white", "red")),
          cluster_columns = FALSE, show_column_names = FALSE, top_annotation = ha2, column_names_gp = gpar(fontsize = 8), show_row_names = FALSE, column_title = "Expression", column_title_gp = gpar(fontsize = 10))+
  Heatmap(FoldChange, name = "Fold Change", col = colorRamp2(c(23.34, 10, 5, 1, 0.31), c("red","purple", "green","blue", "yellow")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(gene_type, name = "Probe Relat.", column_names_gp = gpar(fontsize = 8))+
  Heatmap(anno_gene, name = "Probe Annot.", col = c("TSS200" = "red", "TSS1500" = "blue"), column_names_gp = gpar(fontsize = 8))
 #Heatmap(dist, name = "Distance TSS", col = colorRamp2(c(0, 500, 2000, 15000), c("black", "red", "green", "blue")), column_names_gp = gpar(fontsize = 8))

#pdf("genomic_regions_cluster4_ward2.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive correspondence between methylation, expression and other genomic features", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", legend_title_gp = gpar(fontsize = 8, fontface = "bold"), legend_label_gp = gpar(fontsize = 8))
#dev.off()

#####################################################################
############## ComplexHeatmap with row name in clusters #############
ht_list = Heatmap(mat_meth.impute.TSS.tumor, name = "Methylation", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
                  +                   cluster_columns = column_tree, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 4, column_title = "Methylation", column_title_gp = gpar(fontsize = 10), show_row_names = TRUE, row_names_max_width = unit(4, "cm"), row_title_gp = gpar(fontsize = 12))+
  +   Heatmap(direction, name = "Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  +   Heatmap(mat_expr.impute.TSS.tumor[, column_tree$order], name = "Expression", col = colorRamp2(c(-2.01, 0, 16.87), c("blue","white", "red")),
              +           cluster_columns = FALSE, show_column_names = FALSE, top_annotation = ha2, column_names_gp = gpar(fontsize = 8), show_row_names = FALSE, column_title = "Expression", column_title_gp = gpar(fontsize = 10))+
  +   Heatmap(FoldChange, name = "Fold Change", col = colorRamp2(c(23.34, 10, 5, 1, 0.31), c("red","purple", "green","blue", "yellow")), column_names_gp = gpar(fontsize = 8))+
  +   Heatmap(gene_type, name = "Probe Relat.", column_names_gp = gpar(fontsize = 8))+
  +   Heatmap(anno_gene, name = "Probe Annot.", col = c("TSS200" = "red", "TSS1500" = "blue"), column_names_gp = gpar(fontsize = 8))+
  + Heatmap(dist, name = "Distance TSS", col = colorRamp2(c(0, 500, 2000, 15000), c("black", "red", "green", "blue")), column_names_gp = gpar(fontsize = 8))
> 
  > pdf("genomic_regions_cluster4_ward2.pdf", width = 60, height = 60)
> draw(ht_list, newpage = FALSE, column_title = "Comprehensive correspondence between methylation, expression and other genomic features", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", legend_title_gp = gpar(fontsize = 8, fontface = "bold"), legend_label_gp = gpar(fontsize = 8))
> dev.off()

#####################################################################
############ RNASeq heatmap #########################################
sampleIDs <- colnames(expr.knn)
sampleIDs <- gsub("\\.", "-", sampleIDs)
colnames(expr.knn) <- sampleIDs
samplesDat <- data.frame(matrix(nrow = length(sampleIDs), ncol = 4))
rownames(samplesDat) <- sampleIDs
for (j in 1:length(sampleIDs)) {
  tmpRow <- unlist(strsplit(sampleIDs[j], split = "-"))
  samplesDat[sampleIDs[j], ] <- tmpRow
}
sampleIDs1 <- as.character(samplesDat[, 4])
sampleIDs1 <- substr(sampleIDs1, 1, nchar(sampleIDs1)) 
sampleIDs1 <- as.numeric(sampleIDs1)
normalSamples <- rownames(samplesDat)[sampleIDs1 < 14 & sampleIDs1 > 9]
tumorSamples <- rownames(samplesDat)[sampleIDs1 < 6]
expr.knn <- expr.knn.full[, c(tumorSamples, normalSamples)]


ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)), rep("Normal", length(normalSamples)))), col = list(type = c("Tumor" = "red", "Normal"="blue")))
#ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", length(tumorSamples)))), col = list(type = c("Tumor" = "red")))
column_tree = hclust(dist(t(expr.knn[,1:length(tumorSamples)])), method = "ward.D2")
ht_list = Heatmap(expr.knn[, column_tree$order], name = "Expression", col = colorRamp2(c(-5.25, -1, 0, 3, 5, 13), c("green","blue", "grey", "pink","purple", "red")), 
                  cluster_columns = TRUE, show_column_names = FALSE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 4, show_row_names = FALSE, column_title = "Expression", column_title_gp = gpar(fontsize = 10)) +
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(FoldChange, name = "Fold Change", col = colorRamp2(c(11.5, 0.2, 0.001), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))
pdf("gene_expression_clust4.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential expression analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", legend_title_gp = gpar(fontsize = 8, fontface = "bold"), legend_label_gp = gpar(fontsize = 8))
dev.off()



#####################################################################
# library(ComplexHeatmap)
# library(circlize)
# library(RColorBrewer)
# 
# rand_meth = function(k, mean) {
#   (runif(k) - 0.5)*min(c(1-mean), mean) + mean
# }
# 
# mean_meth = c(rand_meth(300, 0.3), rand_meth(700, 0.7))
# mat_meth = as.data.frame(lapply(mean_meth, function(m) {
#   if(m < 0.3) {
#     c(rand_meth(10, m), rand_meth(10, m + 0.2))
#   } else if(m > 0.7) {
#     c(rand_meth(10, m), rand_meth(10, m - 0.2))
#   } else {
#     c(rand_meth(10, m), rand_meth(10, m + sample(c(1, -1), 1)*0.2))
#   }
#   
# }))
# mat_meth = t(mat_meth)
# rownames(mat_meth) = NULL
# colnames(mat_meth) = paste0("sample", 1:20)
# 
# direction = rowMeans(mat_meth[, 1:10]) - rowMeans(mat_meth[, 11:20])
# direction = ifelse(direction > 0, "hyper", "hypo")
# 
# mat_expr = t(apply(mat_meth, 1, function(x) {
#   x = x + rnorm(length(x), sd = (runif(1)-0.5)*0.4 + 0.5)
#   scale(x)
# }))
# dimnames(mat_expr) = dimnames(mat_meth)
# 
# cor_pvalue = -log10(sapply(seq_len(nrow(mat_meth)), function(i) {
#   cor.test(mat_meth[i, ], mat_expr[i, ])$p.value
# }))
# 
# gene_type = sample(c("protein_coding", "lincRNA", "microRNA", "psedo-gene", "others"), nrow(mat_meth), replace = TRUE, prob = c(6, 1, 1, 1, 1))
# 
# anno_gene = sapply(mean_meth, function(m) {
#   if(m > 0.6) {
#     if(runif(1) < 0.8) return("intragenic")
#   }
#   if(m < 0.3) {
#     if(runif(1) < 0.4) return("TSS")
#   }
#   return("intergenic")
# })
# 
# dist = sapply(mean_meth, function(m) {
#   if(m < 0.6) {
#     if(runif(1) < 0.8) return(round( (runif(1)-0.5)*1000000 + 500000 ))
#   }
#   if(m < 0.3) {
#     if(runif(1) < 0.4) return(round( (runif(1) - 0.5)*1000 + 500))
#   }
#   return(round( (runif(1) - 0.5)*100000 + 50000))
# })
# 
# 
# anno_enhancer_1 = sapply(mean_meth, function(m) {
#   if(m < 0.3) {
#     if(runif(1) < 0.6) return(runif(1))
#   } else if (runif(1) < 0.1) {
#     return(runif(1))
#   } 
#   return(0)
# })
# 
# 
# anno_enhancer_2 = sapply(mean_meth, function(m) {
#   if(m < 0.3) {
#     if(runif(1) < 0.6) return(runif(1))
#   } else if (runif(1) < 0.1) {
#     return(runif(1))
#   } 
#   return(0)
# })
# 
# anno_enhancer_3 = sapply(mean_meth, function(m) {
#   if(m < 0.3) {
#     if(runif(1) < 0.6) return(runif(1))
#   } else if (runif(1) < 0.1) {
#     return(runif(1))
#   } 
#   return(0)
# })
# 
# anno_enhancer = data.frame(enhancer_1 = anno_enhancer_1, enhancer_2 = anno_enhancer_2, enhancer_3 = anno_enhancer_3)
# 
# ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", 10), rep("Control", 10))), col = list(type = c("Tumor" = "red", "Control" = "blue")))
# ha2 = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", 10), rep("Control", 10))), col = list(type = c("Tumor" = "red", "Control" = "blue")), show_legend = FALSE)
# 
# column_tree = hclust(dist(t(mat_meth)))
# 
# ht_list = Heatmap(mat_meth, name = "methylation", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
#                   cluster_columns = column_tree, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 5, column_title = "Methylation", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10)) +
#   Heatmap(direction, name = "direction", col = c("hyper" = "red", "hypo" = "blue"), column_names_gp = gpar(fontsize = 8)) +
#   Heatmap(mat_expr[, column_tree$order], name = "expression", col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
#           cluster_columns = FALSE, top_annotation = ha2, column_names_gp = gpar(fontsize = 8), column_title = "Expression", column_title_gp = gpar(fontsize = 10)) +
#   Heatmap(cor_pvalue, name = "-log10(cor_p)", col = colorRamp2(c(0, 2, 4), c("white", "white", "red")), column_names_gp = gpar(fontsize = 8)) +
#   Heatmap(gene_type, name = "gene type", col = brewer.pal(length(unique(gene_type)), "Set1"), column_names_gp = gpar(fontsize = 8)) +
#   Heatmap(anno_gene, name = "anno_gene", col = brewer.pal(length(unique(anno_gene)), "Set2"), column_names_gp = gpar(fontsize = 8)) +
#   Heatmap(dist, name = "dist_tss", col = colorRamp2(c(0, 10000), c("black", "white")), column_names_gp = gpar(fontsize = 8)) +
#   Heatmap(anno_enhancer, name = "anno_enhancer", col = colorRamp2(c(0, 1), c("white", "orange")), cluster_columns = FALSE, column_names_gp = gpar(fontsize = 8), column_title = "Enhancer", column_title_gp = gpar(fontsize = 10))
# 
# #pdf("genomic_regions.pdf", width = 10, height = 10)
# draw(ht_list, newpage = FALSE, column_title = "Comprehensive correspondence between methylation, expression and other genomic features", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", legend_title_gp = gpar(fontsize = 8, fontface = "bold"), legend_label_gp = gpar(fontsize = 8))
# #dev.off()

