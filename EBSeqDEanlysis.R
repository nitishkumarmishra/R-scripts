library("EBSeq")
#mRNAseq <- read.table(file="PAAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", header=TRUE, sep="\t", row.names=1)
mRNAseq <- read.table(file="mRNAseq_raw_counts.txt", header=TRUE, sep="\t", row.names=1)
#mRNAseq <- mRNAseq[-1,]
sampleIDs <- colnames(mRNAseq)
sampleIDs <- gsub("\\.", "-", sampleIDs)
colnames(mRNAseq) <- sampleIDs
mRNAseq <- data.matrix(mRNAseq)
#samplesDat <- data.frame(matrix(nrow = length(sampleIDs), ncol = 7))
samplesDat <- data.frame(matrix(nrow = length(sampleIDs), ncol = 4))
rownames(samplesDat) <- sampleIDs
for (j in 1:length(sampleIDs)) {
  tmpRow <- unlist(strsplit(sampleIDs[j], split = "-"))
  samplesDat[sampleIDs[j], ] <- tmpRow
}
sampleIDs1 <- as.character(samplesDat[, 4])
#sampleIDs1 <- substr(sampleIDs1, 1, nchar(sampleIDs1) - 1) 
sampleIDs1 <- as.numeric(sampleIDs1)
normalSamples <- rownames(samplesDat)[sampleIDs1 < 20 & sampleIDs1 > 9]
tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]
EBSeqMat <- mRNAseq[, c(normalSamples, tumorSamples)]
Sizes <- MedianNorm(EBSeqMat)
EBOut <- EBTest(Data = EBSeqMat, Conditions = as.factor(c(rep("C1", length(normalSamples)), rep("C2", length(tumorSamples)))), sizeFactors = Sizes, maxround = 5)
EBDERes <- GetDEResults(EBOut, FDR = 0.05)
EBDERes$DEfound## List of DE genes
PlotPostVsRawFC(EBOut, GeneFC)## Plot PosteriorFC Vs. FC
QQP(EBOut) ## QQ plot for checking fitting of model for both C1=Normal and C2=Tumor
DenNHist(EBOut) ## Density plot of q's for empirical vs. simulated
