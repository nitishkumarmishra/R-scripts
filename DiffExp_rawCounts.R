library(DESeq2)
geneExp <- read.csv("CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt", header = TRUE, sep = "\t", row.names = 1)
geneExp <- geneExp[,which(geneExp[1,]=="raw_count")]
geneExp <- geneExp[-c(1),]
colnames(geneExp) <- substr(gsub("\\.", "-", colnames(geneExp)), 1, 16)
cancerID <- grep("01A", colnames(geneExp))
normalID <- grep("11A", colnames(geneExp))
cnts <- cbind(geneExp[,cancerID], geneExp[,normalID])
cnts <- data.matrix(cnts)
#cnts <- cnts[rowSums(cnts) > 0, ]
### Remove sample which have more than 20% zero, similarly sample which have more than 20% zero
cnts <- cnts[apply(cnts,1,function(x) sum(x==0))<ncol(cnts)*0.8,]
cnts <- cnts[apply(cnts,2,function(x) sum(x==0))<nrow(cnts)*0.8,]
cnts <- round(cnts)
cond <- factor(c(rep(2, length(cancerID)), rep(1, length(normalID))))
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
dds <- DESeq(dds)
res <- results(dds, lfcThreshold = 1, pAdjustMethod = "BH")
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj <= 0.01)
tmp <- as.data.frame(resSig)
tmp$geneID <- sapply(strsplit(rownames(tmp),"\\|"),'[[',2)
rownames(tmp) <- sapply(strsplit(rownames(tmp),"\\|"),'[[',1)
tmp$foldchange <- 2^tmp$log2FoldChange

########## voom+limma for diff Exp analysis ########

adj.method = "BH"; adj.pval = 0.01; raw.pval = 0.01; logFC = 1; hmTopUpN = 20; hmTopDownN = 20; meanFilter = 10 
meanCounts <- apply(cnts, 1, mean)
mRNAseq <- data.frame(cnts)
voomMat <- mRNAseq[meanCounts >  meanFilter,]
design <- model.matrix(~0 + factor(c(rep(2, length(cancerID)), rep(1, length(normalID)))))
colnames(design) <- c("Normal", "Tumor")

v <- voom(voomMat, design, plot = TRUE)
fit <- lmFit(v, design)
cont.matrix <- makeContrasts("Tumor-Normal", levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aradeger <- topTable(fit2, adjust.method = adj.method, genelist = fit$genes, number = length(fit2))
aradeger <- data.frame(aradeger[aradeger$adj.P.Val < adj.pval & aradeger$P.Value < raw.pval, ])
#rownames(aradeger) <- aradeger$ID
#aradeger <- aradeger[,-1]
####################################################
v1 <- v  #limma camera for enrichment analysis
geneID <- sapply(strsplit(rownames(aradeger),"\\|"),'[[',2)
idx <- ids2indices(Hs.c2, id=geneID)
rownames(v1) <- sapply(strsplit(rownames(v1),"\\|"),'[[',2)
cam.limma <- camera(v1, idx, design, contrast = cont.matrix, inter.gene.cor = 0.01)
head(cam.limma, 5)
dt <- decideTests(fit2, lfc = 1, adjust.method = "BH", p.value = 0.01)
plotMD(fit2, column=1,status=dt[,1],main=colnames(fit2)[1], xlim=c(-8,13))
#library(Glimma)
#glMDPlot(fit2, coef = 1, status = dt[,1], main=colnames(fit2), ccounts = cnts, samples = colnames(cnts), anno = rownames(cnts))
library(clusterProfiler)
c5 <- read.gmt("C:/Users/nitish.mishra/Desktop/MSigDB/c2.all.v5.1.entrez.gmt")
egmt <- enricher(geneID, TERM2GENE = c5)
head(summary(egmt))
## GSEA from cluster profiler
geneID <- sort(geneID, decreasing = TRUE)
kk2 <- gseKEGG(geneList = geneID, organism = 'hsa', nPerm = 1000, minGSSize = 2, pvalueCutoff = 0.5)
kk <- enrichKEGG(geneID, organism = 'hsa', pvalueCutoff = 0.05)
#######################################################
aradeger <- aradeger[aradeger[, 1] > logFC | aradeger[, 1] < (-1 * logFC), ]
aradeger$geneID <- sapply(strsplit(rownames(aradeger),"\\|"),'[[',2)
rownames(aradeger) <- sapply(strsplit(rownames(aradeger),"\\|"),'[[',1)

volcanoplot(fit2, names = fit2$genes$ID, xlab = "Log Fold Change", ylab = "Log Odds", pch = 16, cex = 0.35)
if (nrow(aradeger) > 2) {
  aradeger <- aradeger[order(aradeger[, 1], decreasing = TRUE), ]
  if (nrow(aradeger) >= (hmTopDownN + hmTopUpN)) {
    if (hmTopUpN > 0) {
      topgenes <- rownames(aradeger)[1:hmTopUpN]
    } else {
      topgenes <- NULL
    }
    if (hmTopDownN > 0) {
      bottomgenes <- rownames(aradeger)[(nrow(aradeger) - (hmTopDownN - 1)):nrow(aradeger)]
    } else {
      bottomgenes <- NULL
    }
    bluered <- colorRampPalette(c("blue", "white", "red"))(256)
    v <- v[c(topgenes, bottomgenes), ]
    v <- apply(v, 2, as.numeric)
    rownames(v) <- c(topgenes, bottomgenes)# colv=NA no clustering of column
    try(heatmap(v, col = bluered, scale = "row", main = "RNASeq"), silent = FALSE)
    #try(heatmap(v, col = bluered, scale = "row", main = "RNASeq", Colv = NA), silent = FALSE)
  } else {
    bluered <- colorRampPalette(c("blue", "white", "red"))(256)
    v <- v[rownames(aradeger), ]
    v <- apply(v, 2, as.numeric)
    rownames(v) <- rownames(aradeger)
    try(heatmap(v, col = bluered, scale = "row", main = "RNASeq", Colv = NA), silent = FALSE)
  }
}

save.image(file = "DiffExp_rawCounts.RData")
