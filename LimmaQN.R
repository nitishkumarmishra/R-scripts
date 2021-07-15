library(limma)
############################
adj.method = "BH"; adj.pval = 0.05; raw.pval = 0.05; logFC = 2; hmTopUpN = 100; hmTopDownN = 100; meanFilter = 10
###########
mRNAseq <- read.table(file="mRNAseq_raw_counts.txt", header=TRUE, sep="\t", row.names=1)
#mRNAseq <- read.table(file="PAAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", header=TRUE, sep="\t", row.names=1)
mRNAseq <- mRNAseq[-1,]
sampleIDs <- colnames(mRNAseq)
sampleIDs <- gsub("\\.", "-", sampleIDs)
colnames(mRNAseq) <- sampleIDs
#mRNAseq <- data.matrix(mRNAseq) ## Convert data frame in matrix
## We have to check the colname of mRNAseq is separated by "-", in not then use gsub to replace "." with "-"
# I found when I am using read.csv or read.table it replace "-" with "." in colname
#I have to also check the length of colnames if its only four e.g "TCGA-2J-AAB4-01" then use ncol = 4 in matrix if colnames is "TCGA-2J-AAB1-01A-11R-A41B-07" then ncol=7
samplesDat <- data.frame(matrix(nrow = length(sampleIDs), ncol = 4))
#samplesDat <- data.frame(matrix(nrow = length(sampleIDs), ncol = 7))
rownames(samplesDat) <- sampleIDs
for (j in 1:length(sampleIDs)) {
  tmpRow <- unlist(strsplit(sampleIDs[j], split = "-"))
  samplesDat[sampleIDs[j], ] <- tmpRow
}
sampleIDs1 <- as.character(samplesDat[, 4])

sampleIDs1 <- as.numeric(sampleIDs1)
normalSamples <- rownames(samplesDat)[sampleIDs1 < 20 & sampleIDs1 > 9]
tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]

############# voom + limma for DE analysis #########
meanCounts <- apply(mRNAseq, 1, mean)
limmaMat <- mRNAseq[meanCounts >  meanFilter, c(normalSamples, tumorSamples)]
design <- model.matrix(~0 + factor(c(rep(1, length(normalSamples)), rep(2, length(tumorSamples)))))
colnames(design) <- c("Normal", "Tumor")
# According to Genome Biology 2013, 14:R95
#library(edgeR)
#nf <- calcNormFactors(mRNAseq)
#v <- voom(limmaMat, design, plot = TRUE, lib.size=colSums(limmaMat) * nf)
#v <- voom(limmaMat, design, plot = TRUE)
#dge <- DGEList(counts=counts)
#limmaMat=log2(limmaMat+1)
v <- voom(limmaMat,design,plot=TRUE,normalize="quantile")
fit <- lmFit(v, design)
cont.matrix <- makeContrasts("Tumor-Normal", levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aradeger <- topTable(fit2, adjust.method = adj.method, genelist = fit$genes, number = length(fit2))
aradeger <- data.frame(aradeger[aradeger$adj.P.Val < adj.pval & aradeger$P.Value < raw.pval, ])
rownames(aradeger) <- aradeger$ID
#aradeger <- aradeger[,-1]
aradeger <- aradeger[aradeger[, 2] > logFC | aradeger[, 2] < (-1 * logFC), ]
#tmpReturn <- new("DGEResult", Dataset = "RNASeq", Toptable = data.frame(aradeger)) ### This line giving error---
#Error in getClass(Class, where = topenv(parent.frame())) :
#"DGEResult" is not a defined class
#listResults <- c(listResults, tmpReturn)
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
    rownames(v) <- c(topgenes, bottomgenes)
    #try(heatmap(v, col = bluered, scale = "row", main = "RNASeq", Colv = NA), silent = FALSE)
    try(heatmap(v, col = bluered, scale = "row", main = "RNASeq", hclustfun=function(d) hclust(d, method="ward.D2")), silent = FALSE)
  } else {
    bluered <- colorRampPalette(c("blue", "white", "red"))(256)
    v <- v[rownames(aradeger), ]
    v <- apply(v, 2, as.numeric)
    rownames(v) <- rownames(aradeger)
    #try(heatmap(v, col = bluered, scale = "row", main = "RNASeq", Colv = NA), silent = FALSE)
    try(heatmap(v, col = bluered, scale = "row", main = "RNASeq", hclustfun=function(d) hclust(d, method="ward.D2")), silent = FALSE)
  }
}

