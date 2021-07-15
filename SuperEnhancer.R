## R program to analyze super enhancer region differential methylation
# Dowload all human super-enhancer data from dbSUPER (database of super-enhancer)
##http://bioinfo.au.tsinghua.edu.cn/dbsuper/download.php
## dbSUPER data in folder "dbSuperEnhancer-hg19" at /storage/gudalab/nmishra/dbSuperEnhancer-hg19
# cat * >dbSuperEnhancer-all-99
# ## Print uniq super-enhancer
# cut -f 1-3 dbSuperEnhancer-all-99 |sort -u >dbSuperEnhancer-all-99-Uniq
## Print duplicate super-enhancer
## cut -f 1-3 dbSuperEnhancer-all-99 |sort |uniq -d
#cut -f 1-4 dbSuperEnhancer-all-99 |egrep -v "SE_63876|SE_64128|SE_63815|SE_63832|SE_63970|SE_57587|SE_22045|SE_63797|SE_63982|SE_64139|SE_63617|SE_64105|SE_63931|SE_63782|SE_63918|SE_64093|SE_63629" >dbSuperEnhancer-all-99-Uniq-1.txt
library(limma)
library(IlluminaHumanMethylation450kprobe)
library(ChAMP)
data("IlluminaHumanMethylation450kprobe")
####################################
dbSUPER <- read.csv("dbSuperEnhancer-all-99-Uniq-1.txt", header = FALSE, sep = "\t")
colnames(dbSUPER) <- c("chr", "start", "end","dbSUPER-ID")
dbSUPER$chr <- gsub("chr", "", dbSUPER$chr)
rownames(dbSUPER) <- dbSUPER$`dbSUPER-ID`
contains <- function(chr, start, end, table) 
  { 
  idx <- table$chr == chr & table$start >= start & table$end <=end
  table[idx,] 
}
n <- dim(dbSUPER)[1]
my.list<- vector("list", n)
for(i in 1:n)
{
  chr1 <- dbSUPER$chr[i]
  start1 <- dbSUPER$start[i]
  end1 <- dbSUPER$end[i]
  dbSUPERID <- dbSUPER$`dbSUPER-ID`[i]
  tmp <- contains(chr1, start1, end1, IlluminaHumanMethylation450kprobe)[,1:6]
  if(dim(tmp)[1]=="0")
  {
    t <- c("NA", "NA", "NA", "NA", "NA", "NA")
    tmp <- rbind(tmp, t)
    rownames(tmp) <- i
    colnames(tmp) <- c("Probe_ID", "chr","strand","start","end","site")
  }
  tmp$dbSUPER <-dbSUPERID
  my.list[[i]] <- tmp
}
finalList <- do.call(rbind, my.list)
finalList1 <- finalList[grep("cg|ch|rs", rownames(finalList)),]

rm(t, i, n, my.list, chr1, start1, end1, tmp, dbSUPERID, finalList)
#### Aggregate methylation data ####
load("CHOL_meth.rda")
Meth <- impute::impute.knn(Meth.common, k=10, rowmax=0.5, colmax=0.8)
Meth <- Meth$data
commonCpGs <- intersect(unique(finalList1$Probe_ID), rownames(Meth))
Meth.common <- Meth[commonCpGs,]
myNorm <- ChAMP::champ.norm(beta = Meth.common, fromIDAT = F, methValue = "B", norm = "BMIQ", fromFile = FALSE, betaFile = Meth.common, filterXY = FALSE, plotBMIQ = TRUE, resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
Meth.common <- myNorm$beta
finalList <- finalList1[(finalList1$Probe_ID%in%commonCpGs),]

#### Remove region which have less than 3 probes #######
name <- dplyr::group_by(finalList, dbSUPER)
name <- dplyr::summarise(name, n())
name <- name[which(name$`n()`>=3),]
dbSUPERID <- name$dbSUPER

#####################################
aggregateMeth <- function(ID, table){
  idx <- table$dbSUPER==ID
  idx <- which(idx)
  table[idx,]
}
n <- length(dbSUPERID)
my.list<- vector("list", n)
for (i in 1:n) {
  ID <- dbSUPERID[i]
  tmp <- aggregateMeth(ID, finalList)
  meth <- Meth.common[rownames(Meth.common)%in%tmp$Probe_ID,]
  meth <- data.frame(cbind(meth, as.character(tmp$dbSUPER)))
  colnames(meth)[dim(meth)[2]] <- c("dbSUPER")
  colnames(meth) <- gsub("\\.","-", colnames(meth))
  meth <- varhandle::unfactor(meth)
  aggre <- aggregate(x = meth[,1:46], by = list(meth$dbSUPER), FUN = mean)
  my.list[[i]]<- aggre
}
finalList2 <- do.call(rbind, my.list)
rownames(finalList2) <- finalList2$Group.1
finalList2$Group.1 <- NULL

#####################################
##### limma Diff Meth analysis ######
Tumor.aggr <- finalList2[,which(substr(colnames(finalList2),14,15)=="01")]
Normal.aggr <- finalList2[,which(substr(colnames(finalList2),14,15)=="11")]
Aggr.Meth <- cbind(Tumor.aggr, Normal.aggr)
Aggr.Meth.M <- wateRmelon::Beta2M(Aggr.Meth)
c1 <- ncol(Tumor.aggr)
c2 <- ncol(Normal.aggr)
design = cbind(Tumor = c(rep(1,  c1), rep(0, c2)), Normal = c(rep(0, c1), rep(1, c2)))
methFit = lmFit(Aggr.Meth.M, design)
contrast.matrix = makeContrasts(TumorVsNormal = Tumor-Normal, levels = design)
methFit2 = contrasts.fit(methFit, contrast.matrix)
methFit2 = eBayes(methFit2)
methRresults = topTable(methFit2, coef = 1, number = Inf, p.value = 0.01, sort.by = "p", resort.by="logFC", adjust.method = "BH")
#####################################
## Annotation of limma result file ##
t.tumor <- as.data.frame(cbind(rowMeans(Tumor.aggr), rowMeans(Normal.aggr), rowMeans(Tumor.aggr) - rowMeans(Normal.aggr), rowMeans(Tumor.aggr) / rowMeans(Normal.aggr), log2(rowMeans(Tumor.aggr) / rowMeans(Normal.aggr))))
colnames(t.tumor) <- c("Tumor", "Normal", "MeanDiff", "FoldChange", "log2FC")
limma.results <- merge(methRresults, t.tumor, by="row.names")
rownames(limma.results) <- limma.results$Row.names
limma.results$Row.names <- NULL
results <- merge(limma.results, dbSUPER, by="row.names")
results <- results[which(results$adj.P.Val <= 0.01  & abs(results$MeanDiff) >= 0.20 & (results$FoldChange >= 2|results$FoldChange <= 0.5)&(results$chr!="X"|results$chr!="Y")),]
rownames(results) <- results$Row.names
results$Row.names<- NULL
results$`dbSUPER-ID` <- NULL

rm(commonCpGs, probe.features, probeInfoALL.lv, Tumor.aggr, Normal.aggr, my.list, myNorm, i, idx, n, Meth, aggre, ID, tmp, name, c1, c2)
#####################################
results <- results[order(results$FoldChange, decreasing = TRUE),]
write.csv(results, file = "DiffDbSuperBMIQ.csv")
save.image(file = "dbSUPER.RData")
