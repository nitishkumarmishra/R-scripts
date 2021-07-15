library(ComplexHeatmap)
library(circlize)
library(NMF)
library(ConsensusClusterPlus)
########################################
load("MethKNN.rda")
Supplemen <- read.csv("Supplementary Table 1.txt", header = TRUE, row.names = 1, sep = "\t")
Supplemen.0.3 <- Supplemen[abs(Supplemen$DeltaBeta)>=0.3,]
mark.0.3 <- rownames(Supplemen.0.3)
Meth-o.3 <- Meth.knn[rownames(Meth.knn)%in%mark.0.3, ]
Meth.knn.0.3.normal <- Meth.knn.0.3[,which(substr(colnames(Meth.knn.0.3),14,15) == "11")]
Meth.knn.0.3.tumor <- Meth.knn.0.3[,which(substr(colnames(Meth.knn.0.3),14,15) == "01")]
estim.r <- nmf(Meth.knn.0.3.tumor, 2:6,  nrun=500, .options='k')
plot(estim.r)
###############
pdf("NMFclusters.pdf", width = 10, height = 12)
consensusmap(estim.r)
dev.off()
#####Select best K number, then calculate silhouette width ####
## Let assume 3 cluster have best result, by checking cophenetic coefficients
si <- silhouette(nmf(Meth.knn.0.3.tumor, 3, nrun = 500))
summary(si) ## To check silhouette width
plot(si)
plot(si, col = c("red", "green", "blue"))

##############################################
pdf("Silhouette.pdf", width = 10, height = 12)
plot(si, col = c("red4", "green4", "blue4"), main="Silhouette Information",do.col.sort=TRUE)
text(x = 0.7, y = 150, paste("Cluster1"), cex = 1.2, col = c("black"))
text(x = 0.7, y = 50, paste("Cluster2"), cex = 1.2, col = c("black"))
text(x = 0.7, y = 5, paste("Cluster3"), cex = 1.2, col = c("black"))
dev.off()
#### ConsensusCluster ######################
title=tempdir()
results = ConsensusClusterPlus(Meth.knn.0.3.tumor,maxK=6,reps=500,pItem=0.8,pFeature=1, title=title,clusterAlg="hc",distance="pearson",seed=1262118388,plot="png")
icl = calcICL(results,title=title,plot="png")
