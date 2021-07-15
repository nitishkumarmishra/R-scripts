####################################################################
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
tmp <- as.matrix(normalSamples.Expression)
tmp.knn <- impute.knn(tmp, k = 15, rowmax = 0.8, colmax = 0.2)
tmp <- tmp.knn$data
tmp1 <- as.matrix(tumorSamples.Expression)
tmp1.knn <- impute.knn(tmp1, k = 15, rowmax = 0.8, colmax = 0.2)
tmp1 <- tmp1.knn$data
x3 <- rowMedians(tmp, na.rm = TRUE)
PAAD.uncv2.mRNAseq_RSEM <- apply(t(tmp1), 1, function(x) t(x)/x3)
rownames(PAAD.uncv2.mRNAseq_RSEM) <- rownames(tmp)
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
#################################################################
################ Overall CNV & gene expression ##################
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
a1.1 <- a1[, !apply(is.na(a1), 1, all)]#Remove col with all NA
a1.1 <- a1[, !apply(is.na(a1), 2, all)]#Remove row with all NA
a2.1 <- a2[, !apply(is.na(a2), 1, all)]
a2.1 <- a2[, !apply(is.na(a2), 2, all)]
x <- intersect(colnames(a1.1),colnames(a2.1))
x1 <- intersect(rownames(a1.1),rownames(a2.1))
a1.1.common <- a1.1[c(x1),c(x)]
a2.1.common <- a2.1[c(x1),c(x)]
t <- apply(a1.1.common, 1, function(x) length(which(!is.na(x))))## Count number of NA in each row
names.less17 <- names(t[which(t<=17)])## Name of gene which have more than 17 NA (i.e 10% sample)
#Total 5795 gene have less than 17 non-NA value, so we used 10702 X 162 matrix for analysis
#In cor.test if less than 17 non-NA it will give error "not enough finite observations"
a1.1.common.1 = a1.1.common[!row.names(a1.1.common)%in%names.less17,]## Remove row with >3 NA
a2.1.common.1 = a2.1.common[!row.names(a2.1.common)%in%names.less17,]
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
a1.1 <- a1[, !apply(is.na(a1), 1, all)]#Remove col with all NA
a1.1 <- a1[, !apply(is.na(a1), 2, all)]#Remove row with all NA
a2.1 <- a2[, !apply(is.na(a2), 1, all)]
a2.1 <- a2[, !apply(is.na(a2), 2, all)]
x <- intersect(colnames(a1.1),colnames(a2.1))
x1 <- intersect(rownames(a1.1),rownames(a2.1))
a1.1.common <- a1.1[c(x1),c(x)]
a2.1.common <- a2.1[c(x1),c(x)]
t <- apply(a1.1.common, 1, function(x) length(which(!is.na(x))))## Count number of NA in each row
names.less17 <- names(t[which(t<=17)])## Name of gene which have more than 17 NA (i.e 10% sample)
#Total 5795 gene have less than 17 non-NA value, so we used 10702 X 162 matrix for analysis
#In cor.test if less than 17 non-NA it will give error "not enough finite observations"
m1.1.common.1 = a1.1.common[!row.names(a1.1.common)%in%names.less17,]## Remove row with >3 NA
m2.1.common.1 = a2.1.common[!row.names(a2.1.common)%in%names.less17,]
rm(a, i, j, tmp, tmp1, tmpRow, x3, sampleIDs1, samplesDat, sampleIDs, colName, tumorSamples, normalSamples, a1, a1.1, a1.1.common, t, names.less17, x, x1, a2, a2.1, a2.1.common)

################# Spearman's rank correlation ###################
################# For CNV loss ###################
ll <- mapply(function(x,y)cor.test(a1.1.common.1[x,],a2.1.common.1[y,], method = "spearman", alternative = "greater"),
             1:nrow(a1.1.common.1),
             1:nrow(a2.1.common.1),
             SIMPLIFY=FALSE)
# ll <- mapply(function(x,y)cor.test(all_data.by_genes.common[x,],PAAD.uncv2.mRNAseq_RSEM.common.impute[y,], method = "spearman", alternative = "greater"),
#              1:nrow(all_data.by_genes.common),
#              1:nrow(PAAD.uncv2.mRNAseq_RSEM.common),
#              SIMPLIFY=FALSE)
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
CNVcommon.p <- cbind(p.value, cor.value)
rownames(CNVcommon.p) <- rownames(a1.1.common.1)
CNVcommon.p <- data.frame(CNVcommon.p)
adjust <- mt.rawp2adjp(CNVcommon.p$p.value, proc=(c("Bonferroni")))
adjust.order <- adjust$adjp[order(adjust$index),]
adjust.order.1 <- data.frame(adjust.order)
final.CNVcommon <- cbind(CNVcommon.p, adjust.order.1)
CNVloss.adjP.0.01 <- final.CNVcommon[(final.CNVcommon$rawp <= 0.01) & (final.CNVcommon$Bonferroni <= 0.01),]
dim(CNVloss.adjP.0.01)
CNVloss.adjP.0.01[order(rownames(CNVloss.adjP.0.01), decreasing = FALSE),]### Order by Gene name
CNVloss.adjP.0.01[order(CNVloss.adjP.0.01[,4], decreasing = FALSE),]## Order by Bonferroni adjusted p-value
rm(ll, cor.value, p.value, CNVcommon.p, adjust, adjust.order, adjust.order.1)
################# For CNV gain ###################
ll <- mapply(function(x,y)cor.test(m1.1.common.1[x,],m2.1.common.1[y,], method = "spearman", alternative = "greater"),
             1:nrow(m1.1.common.1),
             1:nrow(m2.1.common.1),
             SIMPLIFY=FALSE)
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
CNVcommon.p <- cbind(p.value, cor.value)
rownames(CNVcommon.p) <- rownames(m1.1.common.1)
CNVcommon.p <- data.frame(CNVcommon.p)
adjust <- mt.rawp2adjp(CNVcommon.p$p.value, proc=(c("Bonferroni")))
adjust.order <- adjust$adjp[order(adjust$index),]
adjust.order.1 <- data.frame(adjust.order)
final.CNVcommon <- cbind(CNVcommon.p, adjust.order.1)
CNVgain.adjP.0.01 <- final.CNVcommon[(final.CNVcommon$rawp <= 0.01) & (final.CNVcommon$Bonferroni <= 0.01),]
dim(CNVgain.adjP.0.01)
CNVgain.adjP.0.01[order(rownames(CNVgain.adjP.0.01), decreasing = FALSE),]### Order by Gene name
CNVgain.adjP.0.01[order(CNVgain.adjP.0.01[,4], decreasing = FALSE),]#
rm(ll, cor.value, p.value, CNVcommon.p, adjust, adjust.order, adjust.order.1)
#######################################################

