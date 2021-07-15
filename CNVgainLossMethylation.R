## R program to calculate Spearman's correlation between gene expression and CNV
# No concept of gain & loss, this program common genes and samples in CNV and gene expression
# Calculate Spearman's correlation between CNV and gene expression
library(matrixStats)
library(multtest)
library(impute)
###############################################
all_threshold.by_genes <- read.csv("all_thresholded.by_genes.txt", header = TRUE, row.names = 1, sep = "\t")
sampleIDs <- colnames(all_threshold.by_genes)
sampleIDs <- gsub("\\.", "-", sampleIDs)
sampleIDs <- substring(sampleIDs, 1, 15) #This command will print only "TCGA-2J-AAB1-01" in name of TCGA samples
colnames(all_threshold.by_genes) <- sampleIDs
all_threshold.by_genes$`Locus-ID`<- NULL
all_threshold.by_genes$Cytoband <- NULL
#################################################
all_data.by_genes <- read.csv("all_data_by_genes.txt", header = TRUE, row.names = 1, sep = "\t")
sampleIDs <- colnames(all_data.by_genes)
sampleIDs <- gsub("\\.", "-", sampleIDs)
sampleIDs <- substring(sampleIDs, 1, 15) #This command will print only "TCGA-2J-AAB1-01" in name of TCGA samples
colnames(all_data.by_genes) <- sampleIDs
all_data.by_genes$`Locus-ID`<- NULL
all_data.by_genes$Cytoband <- NULL
############### Load expression data #############
PAAD.uncv2.mRNAseq_RSEM_normalized_log2 <- read.csv("PAAD.uncv2.mRNAseq_RSEM_normalized_log2.txt", header = TRUE, row.names = 1, sep = "\t")
colName <- colnames(PAAD.uncv2.mRNAseq_RSEM_normalized_log2)
colName <- gsub("\\.", "-", colName)
colnames(PAAD.uncv2.mRNAseq_RSEM_normalized_log2) <- colName
PAAD.uncv2.mRNAseq_RSEM_normalized_log2$`gene-1` <- NULL
sampleIDs <- colnames(PAAD.uncv2.mRNAseq_RSEM_normalized_log2)
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
tumorSamples <- rownames(samplesDat)[sampleIDs1 < 3]
normalSamples.Expression <- PAAD.uncv2.mRNAseq_RSEM_normalized_log2[, c(normalSamples)]
tumorSamples.Expression <- PAAD.uncv2.mRNAseq_RSEM_normalized_log2[, c(tumorSamples)]
# tmp <- as.matrix(normalSamples.Expression)
# tmp1 <- as.matrix(tumorSamples.Expression)
# x3 <- rowMedians(tmp, na.rm = TRUE)
tmp <- as.matrix(normalSamples.Expression)
tmp.knn <- impute.knn(tmp, k = 15, rowmax = 0.8, colmax = 0.2)
tmp <- tmp.knn$data
tmp1 <- as.matrix(tumorSamples.Expression)
tmp1.knn <- impute.knn(tmp1, k = 15, rowmax = 0.8, colmax = 0.2)
tmp1 <- tmp1.knn$data
x3 <- rowMedians(tmp, na.rm = TRUE)
PAAD.uncv2.mRNAseq_RSEM <- apply(t(tmp1), 1, function(x) t(x)/x3)
rownames(PAAD.uncv2.mRNAseq_RSEM) <- rownames(tmp)
### Function for expression data normalization with normal samples ######
# The commented lines for prgram without normalization with mean of normal samples
##########################################
x <- intersect(colnames(all_data.by_genes),colnames(PAAD.uncv2.mRNAseq_RSEM_normalized_log2))
x1 <- intersect(rownames(all_data.by_genes), rownames(PAAD.uncv2.mRNAseq_RSEM_normalized_log2))
Expr.common <- PAAD.uncv2.mRNAseq_RSEM[c(x1),c(x)]
# PAAD.uncv2.mRNAseq_RSEM.common <- PAAD.uncv2.mRNAseq_RSEM_normalized_log2[c(x1),c(x)]
# PAAD.uncv2.mRNAseq_RSEM.common <- as.matrix(PAAD.uncv2.mRNAseq_RSEM.common)
all_data.by_genes.common <- all_data.by_genes[c(x1),c(x)]
# PAAD.uncv2.mRNAseq_RSEM.common.impute <- impute.knn(PAAD.uncv2.mRNAseq_RSEM.common, k = 15, rowmax = 0.8, colmax = 0.8)
# Expr.common <- PAAD.uncv2.mRNAseq_RSEM.common.impute$data
all_threshold.common <- all_threshold.by_genes[c(x1),c(x)]
CNV.common <- as.matrix(all_data.by_genes.common)
##########################################
# m <- matrix( sample(12,100,replace = TRUE) , 25 , 4 )
# n <- matrix( sample(12,100,replace = TRUE) , 25 , 4 )
# colnames(m) <- c("one","two","three","four")
#CNV.loss <- apply(all_threshold.common,1, function(x) which(x <= -1))
## Make two list of index
#m1 <- apply(m,1,function(x) which(x >=2))
###########################################
#
#########################################
ll <- suppressWarnings(mapply(function(x,y)cor.test(Expr.common[x,],CNV.common[y,], method = "spearman", alternative = "greater"),
1:nrow(CNV.common),
1:nrow(Expr.common),
SIMPLIFY=FALSE))
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
CNVcommon.p <- cbind(cor.value, p.value)
rownames(CNVcommon.p) <- rownames(CNV.common)
CNVcommon.p <- data.frame(CNVcommon.p)
adjust <- mt.rawp2adjp(CNVcommon.p$p.value, proc=(c("Bonferroni")))
adjust.order <- adjust$adjp[order(adjust$index),]
adjust.order.1 <- data.frame(adjust.order)
final.CNVcommon <- cbind(CNVcommon.p, adjust.order.1)
CNVcommon.adjP.0.01 <- final.CNVcommon[(final.CNVcommon$rawp <= 0.01) & (final.CNVcommon$Bonferroni <= 0.01),]
dim(CNVcommon.adjP.0.01)
CNVcommon.adjP.0.01[order(CNVcommon.adjP.0.01[,4], decreasing = FALSE),]
#######################################################
rm(x, x1, j, ll, tmp, tmp1, tmpRow, x3, sampleIDs1, samplesDat, sampleIDs, colName, tumorSamples, normalSamples, cor.value, p.value, adjust, adjust.order, final.CNVcommon)
