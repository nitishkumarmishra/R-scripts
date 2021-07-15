#### Pipeline for miRNA-mRNA integration by using mirComb tool
## In this pipeline, we used log2 of reads and limma for DE analysis
library(miRComb)
miRNAExp <- miRNAExp[rowSums(miRNAExp) >= 10, ]
geneExp <- geneExp[rowSums(geneExp) > 10, ]
cancerID <- grep("01A", colnames(miRNAExp))
normalID <- grep("11A", colnames(miRNAExp))
miRNAExp <- cbind(miRNAExp[,normalID], miRNAExp[,cancerID])
cancerID <- grep("01A", colnames(geneExp))
normalID <- grep("11A", colnames(geneExp))
geneExp <- cbind(geneExp[,normalID], geneExp[,cancerID])
####################################################################
### log2 transformation of data
mRNA <- log2(geneExp+1)
miRNA <- log2(miRNAExp+1)
pheno.mRNA <- data.frame(group = c(rep("H", length(normalID)), rep("D", length(cancerID))),
                         DvH = as.numeric(c(rep("0", length(normalID)), rep("1", length(cancerID)))))
pheno.miRNA <- data.frame(group = c(rep("H", length(normalID)), rep("D", length(cancerID))),
                          DvH = as.numeric(c(rep("0", length(normalID)), rep("1", length(cancerID)))))
rownames(pheno.mRNA) <- colnames(mRNA)
rownames(pheno.miRNA) <- colnames(miRNA)
mRNA <- as.matrix(mRNA)
miRNA <- as.matrix(miRNA)
#####################################################################
data.obj<-new("corObject",dat.miRNA=as.matrix(miRNA),dat.mRNA=as.matrix(mRNA), pheno.miRNA=pheno.miRNA,pheno.mRNA=pheno.mRNA)
### Plot to check data
plotCordist(data.obj,subset="mRNA",type="dist")
plotCordist(data.obj,subset="miRNA",type="dist")
boxplotSamples(data.obj,subset="mRNA")
boxplotSamples(data.obj,subset="miRNA")
plotHeatmap(data.obj,"mRNA")
plotHeatmap(data.obj,"miRNA")
plotPca(data.obj,subset="mRNA")
plotPca(data.obj,subset="miRNA")
#####################################################################
data.obj<-addDiffexp(data.obj,"miRNA",classes="DvH",method.dif="limma")
data.obj<-addDiffexp(data.obj,"mRNA",classes="DvH",method.dif="limma")
data.obj<-addSig(data.obj,"mRNA",adj.pval=0.05,FC=1.5)
data.obj<-addSig(data.obj,"miRNA",adj.pval=0.05)
data.obj<-addCorrelation(data.obj,alternative="less")
data.obj<-addNet(data.obj)
data(microCosm_v5_18)
data(targetScan_v6.2_18)
data.obj<-addDatabase(data.obj,database=c("microCosm_v5_18","targetScan_v6.2_18"))
data.obj<-correctPval(data.obj, pval="pval",method.adj="BH")
data.obj<-addScore(data.obj)
plotCorrelation(data.obj,miRNA="hsa-let-7i",mRNA="A1CF",type="cor",col.color="group",sample.names=TRUE)
plotCorrelation(data.obj,miRNA="hsa-let-7i",mRNA="A1CF",type="residuals")