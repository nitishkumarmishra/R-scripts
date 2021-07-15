library("edgeR")
library("impute")
mRNAseq1 <- read.table(file="mRNAseq_raw_counts.txt", header=TRUE, sep="\t", row.names=1)
#mRNAseq <- mRNAseq[-1,]
mRNAseq1<-as.matrix(mRNAseq1)
mRNAseq.knn <- impute.knn(mRNAseq1, k=15, rowmax = 0.2, colmax = 0.2)
mRNAseq <- mRNAseq.knn$data
sampleIDs <- colnames(mRNAseq)
sampleIDs <- gsub("\\.", "-", sampleIDs)
colnames(mRNAseq) <- sampleIDs
mRNAseq <- data.matrix(mRNAseq)
#mRNAseq <- round(mRNAseq)
samplesDat <- data.frame(matrix(nrow = length(sampleIDs), ncol = 4))
rownames(samplesDat) <- sampleIDs
for (j in 1:length(sampleIDs)) {
  tmpRow <- unlist(strsplit(sampleIDs[j], split = "-"))
  samplesDat[sampleIDs[j], ] <- tmpRow
}
sampleIDs1 <- as.character(samplesDat[, 4])
#sampleIDs1 <- substr(sampleIDs1, 1, nchar(sampleIDs1) - 1)
sampleIDs1 <- as.numeric(sampleIDs1)
normalSamples <- rownames(samplesDat)[sampleIDs1 < 15 & sampleIDs1 > 9]
#tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]
tumorSamples <- rownames(samplesDat)[sampleIDs1 < 6]
mRNAseq <- mRNAseq[, c(normalSamples, tumorSamples)]
mRNAseq <- round(mRNAseq, 0)
edgeRMat<- DGEList(counts=mRNAseq, group=as.factor(c(rep("1", length(normalSamples)), rep("2", length(tumorSamples)))))
edgeRMat<- calcNormFactors(edgeRMat)
edgeRMat <- estimateCommonDisp(edgeRMat)
edgeRMat <- estimateTagwiseDisp(edgeRMat)
et <- exactTest(edgeRMat)
res <- topTags(et,n=nrow(et))
aradeger <- res$table[res$table[, 1] > 2.0 | res$table[, 1] < (-1 * 2.0), ]
#aradeger <- res$table[res$table[, 1] > 1.5 | res$table[, 1] < (-1 * 1.5), ]
aradeger1 <- aradeger[(aradeger[,3] <= 0.05) & (aradeger[,4] <= 0.05),]
#aradeger1 <- aradeger[aradeger[,3] <= 0.05,]
#aradeger2 <- aradeger1[aradeger1[,4] <= 0.05,]
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))
detags <- rownames(edgeRMat)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h = c(-2, 2), col = "blue")