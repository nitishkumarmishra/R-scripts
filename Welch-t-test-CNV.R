## Welch's test for gene expression correlation between two group
## Group Set1- CNVgain Set2- CNVloss set3- CNVunchanged
library(matrixStats)
library(multtest)
library(impute)
############## Load CNV data value ###############
all_data.by_genes <- read.csv("all_data_by_genes.txt", header = TRUE, row.names = 1, sep = "\t")
sampleIDs <- colnames(all_data.by_genes)
sampleIDs <- gsub("\\.", "-", sampleIDs)
sampleIDs <- substring(sampleIDs, 1, 15) #This command will print only "TCGA-2J-AAB1-01" in name of TCGA samples
colnames(all_data.by_genes) <- sampleIDs
all_data.by_genes$`Locus-ID`<- NULL
all_data.by_genes$Cytoband <- NULL
############ Load CNV data threshold #############
all_thresholded.by_genes <- read.csv("all_thresholded.by_genes.txt", header = TRUE, row.names = 1, sep = "\t")
sampleIDs <- colnames(all_thresholded.by_genes)
sampleIDs <- gsub("\\.", "-", sampleIDs)
sampleIDs <- substring(sampleIDs, 1, 15) #This command will print only "TCGA-2J-AAB1-01" in name of TCGA samples
colnames(all_thresholded.by_genes) <- sampleIDs
all_thresholded.by_genes$`Locus-ID`<- NULL
all_thresholded.by_genes$Cytoband <- NULL
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
#tmp <- as.matrix(normalSamples.Expression)
#tmp.knn <- impute.knn(tmp, k = 15, rowmax = 0.8, colmax = 0.2)
#tmp <- tmp.knn$data
tmp1 <- as.matrix(tumorSamples.Expression)
tmp1.knn <- impute.knn(tmp1, k = 15, rowmax = 0.8, colmax = 0.2)
tmp1 <- tmp1.knn$data
#x3 <- rowMedians(tmp, na.rm = TRUE)
#PAAD.uncv2.mRNAseq_RSEM <- apply(t(tmp1), 1, function(x) t(x)/x3)
PAAD.uncv2.mRNAseq_RSEM <- tmp1
rownames(PAAD.uncv2.mRNAseq_RSEM) <- rownames(tmp1)
rm(tmp, tmp1, x3, tmp1.knn, tmp.knn)
############### Matrix for common gene in CNV and expression ###############
x <- intersect(colnames(all_data.by_genes),colnames(PAAD.uncv2.mRNAseq_RSEM))
x1 <- intersect(rownames(all_data.by_genes), rownames(PAAD.uncv2.mRNAseq_RSEM))
Expr.common <- PAAD.uncv2.mRNAseq_RSEM[c(x1),c(x)]
all_data.by_genes.common <- all_data.by_genes[c(x1),c(x)]
all_data.by_genes.common <- all_data.by_genes[c(x1),c(x)]
all_thresholded.by_genes.common <- all_thresholded.by_genes[c(x1),c(x)]
CNV.common <- as.matrix(all_data.by_genes.common)
CNV.loss <- apply(all_thresholded.by_genes.common,1, function(x) which(x <= -1.0))
CNV.gain <- apply(all_thresholded.by_genes.common,1, function(x) which(x >= 1.0))
CNV.unchanged <- apply(all_thresholded.by_genes.common,1, function(x) which(x == 0))
#################################################################
################# File for the CNV unchanged ####################
a1 <- matrix(nrow = nrow(CNV.common), ncol = ncol(CNV.common))
a2 <- matrix(nrow = nrow(CNV.common), ncol = ncol(CNV.common))
rownames(a1) <- rownames(CNV.common)
rownames(a2) <- rownames(CNV.common)
colnames(a1) <- colnames(CNV.common)
colnames(a2) <- colnames(CNV.common)
for(i in 1:nrow(CNV.common))
{
  a <- unlist(CNV.unchanged[i])
  a1[i,a] <- CNV.common[i,a]
  a2[i,a] <- Expr.common[i,a]
}
#a1.1 <- a1[, !apply(is.na(a1), 1, all)]#Remove col with all NA
a1.1 <- a1[, !apply(is.na(a1), 2, all)]#Remove row with all NA
#a2.1 <- a2[, !apply(is.na(a2), 1, all)]
a2.1 <- a2[, !apply(is.na(a2), 2, all)]
x <- intersect(colnames(a1.1),colnames(a2.1))
x1 <- intersect(rownames(a1.1),rownames(a2.1))
a1.1.common <- a1.1[c(x1),c(x)]
a2.1.common <- a2.1[c(x1),c(x)]
#t <- apply(a2.1.common, 1, function(x) length(which(!is.na(x))))## Count number of NA in each row
#names.less17 <- names(t[which(t<=17)])## Name of gene which have more than 17 NA (i.e 10% sample)
#Total 5795 gene have less than 17 non-NA value, so we used 10702 X 162 matrix for analysis
#In cor.test if less than 17 non-NA it will give error "not enough finite observations"
#a1.1.common.1 = a1.1.common[!row.names(a1.1.common)%in%names.less17,]## Remove row with >3 NA
#a2.1.common.1 = a2.1.common[!row.names(a2.1.common)%in%names.less17,]
Expr.common.unchanged <- a2.1.common
############### File for the CNV loss analysis ##################
a1 <- matrix(nrow = nrow(CNV.common), ncol = ncol(CNV.common))
a2 <- matrix(nrow = nrow(CNV.common), ncol = ncol(CNV.common))
rownames(a1) <- rownames(CNV.common)
rownames(a2) <- rownames(CNV.common)
colnames(a1) <- colnames(CNV.common)
colnames(a2) <- colnames(CNV.common)
for(i in 1:nrow(CNV.common))
{
  a <- unlist(CNV.loss[i])
  a1[i,a] <- CNV.common[i,a]
  a2[i,a] <- Expr.common[i,a]
}
#a1.1 <- a1[, !apply(is.na(a1), 1, all)]#Remove col with all NA
a1.1 <- a1[, !apply(is.na(a1), 2, all)]#Remove row with all NA
#a2.1 <- a2[, !apply(is.na(a2), 1, all)]
a2.1 <- a2[, !apply(is.na(a2), 2, all)]
x <- intersect(colnames(a2.1),colnames(a2.1))
x1 <- intersect(rownames(a2.1),rownames(a2.1))
a1.1.common <- a1.1[c(x1),c(x)]
a2.1.common <- a2.1[c(x1),c(x)]
#t <- apply(a2.1.common, 1, function(x) length(which(!is.na(x))))## Count number of NA in each row
#names.less17 <- names(t[which(t<=17)])## Name of gene which have more than 17 NA (i.e 10% sample)
#Total 5795 gene have less than 17 non-NA value, so we used 10702 X 162 matrix for analysis
#In cor.test if less than 17 non-NA it will give error "not enough finite observations"
#a1.1.common.1 = a1.1.common[!row.names(a1.1.common)%in%names.less17,]## Remove row with >3 NA
#a2.1.common.1 = a2.1.common[!row.names(a2.1.common)%in%names.less17,]
Expr.common.loss <- a2.1.common
rm(a, i, j, tmp, tmp1, tmpRow, x3, sampleIDs1, samplesDat, sampleIDs, colName, tumorSamples, normalSamples, a1, a1.1, a1.1.common, t, names.less17, x, x1, a2, a2.1, a2.1.common)
#################################################################
################ File for the CNV gain analysis #################
a1 <- matrix(nrow = nrow(CNV.common), ncol = ncol(CNV.common))
a2 <- matrix(nrow = nrow(CNV.common), ncol = ncol(CNV.common))
rownames(a1) <- rownames(CNV.common)
rownames(a2) <- rownames(CNV.common)
colnames(a1) <- colnames(CNV.common)
colnames(a2) <- colnames(CNV.common)
for(i in 1:nrow(CNV.common))
{
  a <- unlist(CNV.gain[i])
  a1[i,a] <- CNV.common[i,a]
  a2[i,a] <- Expr.common[i,a]
}
#a1.1 <- a1[, !apply(is.na(a1), 1, all)]#Remove col with all NA
a1.1 <- a1[, !apply(is.na(a1), 2, all)]#Remove row with all NA
#a2.1 <- a2[, !apply(is.na(a2), 1, all)]
a2.1 <- a2[, !apply(is.na(a2), 2, all)]
x <- intersect(colnames(a1.1),colnames(a2.1))
x1 <- intersect(rownames(a1.1),rownames(a2.1))
a1.1.common <- a1.1[c(x1),c(x)]
a2.1.common <- a2.1[c(x1),c(x)]
#t <- apply(a2.1.common, 1, function(x) length(which(!is.na(x))))## Count number of NA in each row
#names.less17 <- names(t[which(t<=17)])## Name of gene which have more than 17 NA (i.e 10% sample)
#Total 5795 gene have less than 17 non-NA value, so we used 10702 X 162 matrix for analysis
#In cor.test if less than 17 non-NA it will give error "not enough finite observations"
#m1.1.common.1 = a1.1.common[!row.names(a1.1.common)%in%names.less17,]## Remove row with >3 NA
#m2.1.common.1 = a2.1.common[!row.names(a2.1.common)%in%names.less17,]
Expr.common.gain <- a2.1.common
rm(a, i, j, tmp, tmp1, tmpRow, x3, sampleIDs1, samplesDat, sampleIDs, colName, tumorSamples, normalSamples, a1, a1.1, a1.1.common, t, names.less17, x, x1, a2, a2.1, a2.1.common)

################# Welch's t-test ###################
################# For CNV gain ###################
ll <- mapply(function(x,y)t.test(Expr.common.gain[x,],Expr.common.unchanged[y,], alternative = "g"),
             1:nrow(Expr.common.gain),
             1:nrow(Expr.common.unchanged),
             SIMPLIFY=FALSE)
t.statistic <- sapply(ll,'[[','statistic')
p.value <- sapply(ll,'[[','p.value')
CNVcommon.p <- cbind(p.value, t.statistic)
rownames(CNVcommon.p) <- rownames(Expr.common.gain)
CNVcommon.p <- data.frame(CNVcommon.p)
adjust <- mt.rawp2adjp(CNVcommon.p$p.value, proc=(c("Bonferroni")))
adjust.order <- adjust$adjp[order(adjust$index),]
adjust.order.1 <- data.frame(adjust.order)
final.CNVcommon <- cbind(CNVcommon.p, adjust.order.1)
Welch.CNVgain.adjP.0.05 <- final.CNVcommon[(final.CNVcommon$rawp <= 0.05) & (final.CNVcommon$Bonferroni <= 0.05),]
################# For CNV loss ###################
ll <- mapply(function(x,y)t.test(Expr.common.loss[x,],Expr.common.unchanged[y,], alternative = "l"),
             1:nrow(Expr.common.gain),
             1:nrow(Expr.common.unchanged),
             SIMPLIFY=FALSE)
t.statistic <- sapply(ll,'[[','statistic')
p.value <- sapply(ll,'[[','p.value')
CNVcommon.p <- cbind(p.value, t.statistic)
rownames(CNVcommon.p) <- rownames(Expr.common.gain)
CNVcommon.p <- data.frame(CNVcommon.p)
adjust <- mt.rawp2adjp(CNVcommon.p$p.value, proc=(c("Bonferroni")))
adjust.order <- adjust$adjp[order(adjust$index),]
adjust.order.1 <- data.frame(adjust.order)
final.CNVcommon <- cbind(CNVcommon.p, adjust.order.1)
Welch.CNVloss.adjP.0.05 <- final.CNVcommon[(final.CNVcommon$rawp <= 0.05) & (final.CNVcommon$Bonferroni <= 0.05),]
##########################
rm(ll, adjust, t.statistic, p.value, CNV.common, adjust.order, adjust.order.1, all_data.by_genes, all_thresholded.by_genes, normalSamples.Expression, tumorSamples.Expression, CNV.gain, CNV.loss, CNV.unchanged, CNV.common, PAAD.uncv2.mRNAseq_RSEM, PAAD.uncv2.mRNAseq_RSEM_normalized_log2, all_data.by_genes.common, all_thresholded.by_genes.common, final.CNVcommon, Expr.common, Expr.common.gain, Expr.common.loss, Expr.common.unchanged, CNVcommon.p)
