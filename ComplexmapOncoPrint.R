library(ComplexHeatmap)
library(circlize)
######### R program to make file for non-silent mutation
####### Here Non-silent Mutect file is file which have Somatic mutation
#De_novo_Start_InFrame, De_novo_Start_OutOfFrame, Silent, 
Mutect.NonSilent <- read.csv("NonSilentSNV", header = TRUE, sep = "\t")
################################
NMF.cluster <- si
rownames(NMF.cluster) <- substr(colnames(Meth.knn.0.3.tumor),1,16)
write.csv(NMF.cluster, file = "NMFcluster.csv")
NMF.cluster <- read.csv("NMFcluster.csv", header = TRUE, row.names = 1)
NMF.cluster.order <- NMF.cluster[order(NMF.cluster$cluster, decreasing = FALSE),]
colnames(Meth.knn.0.3.tumor) <- substr(colnames(Meth.knn.0.3.tumor),1,16)
Meth.knn.0.3.tumor.NMFcluster <- Meth.knn.0.3.tumor[,c(rownames(NMF.cluster.order))]
###############################
TP53 <- Mutect.NonSilent[which(as.character(Mutect.NonSilent$Hugo_Symbol)=="TP53"),]
tumor.sample <- colnames(Meth.knn.0.3.tumor.NMFcluster)
#tumor.sample <- substr(tumor.sample, 1,16)
TP53.sample <- substr(TP53$Tumor_Sample_Barcode, 1,16)
ifelse(tumor.sample%in%TP53.sample,1,0)
################################
# t <- c(1:15)
# m <- median(t)
# sd <- sd(t)
# ifelse(t >= m+sd, "2", ifelse(t < m-sd, "1","0"))
# ifelse(t >= m+sd, "High", ifelse(t < m-sd, "Low","Modr"))

################################
CNV_all_data_by_genes <- read.csv("CNV_all_data_by_genes.txt", header = TRUE, row.names = 1, sep = "\t")
CNV_all_data_by_genes$Locus.ID <- NULL
CNV_all_data_by_genes$Cytoband <- NULL
sampleID <- gsub("\\.","-",colnames(CNV_all_data_by_genes))
colnames(CNV_all_data_by_genes) <- substr(sampleID, 1,16)
TP53.cnv <- CNV_all_data_by_genes[which(rownames(CNV_all_data_by_genes)=="TP53"),]
TP53.amplification <- colnames(TP53.cnv[(which(TP53.cnv >= 1.32))])## (2^1.32)*2 = 4.993
TP53.deletion <- colnames(TP53.cnv[(which(TP53.cnv <= -0.2))])
TP53.final <- ifelse(tumor.sample%in%TP53.sample,"SNV", ifelse(tumor.sample%in%TP53.amplification, "Amplification", ifelse(tumor.sample%in%TP53.deletion,"Deletion","NA")))
## Use NA in place of blank space
#################################
## Here Just change gene name in these lines and generate file for that genes #######
# tumor.sample <- colnames(Meth.knn.0.3.tumor.NMFcluster)
# KMT2C <- Mutect.NonSilent[which(as.character(Mutect.NonSilent$Hugo_Symbol)=="KMT2C"),]
# KMT2C.sample <- substr(KMT2C$Tumor_Sample_Barcode, 1,16)
# KMT2C.cnv <- CNV_all_data_by_genes[which(rownames(CNV_all_data_by_genes)=="KMT2C"),]
# KMT2C.amplification <- colnames(KMT2C.cnv[(which(KMT2C.cnv >= 1.32))])
# KMT2C.deletion <- colnames(KMT2C.cnv[(which(KMT2C.cnv <= -0.2))])
# KMT2C.final <- ifelse(tumor.sample%in%KMT2C.sample,"SNV", ifelse(tumor.sample%in%KMT2C.amplification, "Amplification", ifelse(tumor.sample%in%KMT2C.deletion,"Deletion","NA")))
# ##########################
ha <- HeatmapAnnotation(Clusters = c(rep("Cluster1", 43), rep("Cluster2",111), rep("Cluster3",30)), col = list(Clusters = c("Cluster1" = "blue4", "Cluster2"="green4", "Cluster3"="red3")))
Heatmap(Meth.knn.0.3.tumor.NMFcluster,name = "Beta value", top_annotation = ha,col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), cluster_columns = FALSE, cluster_rows = TRUE, show_row_names = FALSE, km = 5,)
#Heatmap(Meth.knn.0.3.tumor.NMFcluster, col = colorRamp2(c(0, 0.2, 0.4, 0.60, 0.8, 1), c("blue4", "green", "gray", "violetred", "red2","red4")), cluster_columns = FALSE, cluster_rows = TRUE, show_row_names = FALSE, km = 5, )
#draw(ha,newpage = FALSE, column_title = "Clsuering of DNA methylation")
# ###############
# df <- data.frame(Clusters = c(rep("Cluster1", 43), rep("Cluster2",111), rep("Cluster3",30)), TP53=c(TP53.final))
# ha <- HeatmapAnnotation(df = df, col = list(Clusters = c("Cluster1" = "blue4", "Cluster2"="green4", "Cluster3"="red3")), gap = unit(c(2, 4), "mm"))
# Heatmap(Meth.knn.0.3.tumor.NMFcluster,name = "Beta value", top_annotation = ha,col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), cluster_columns = FALSE, cluster_rows = TRUE, show_row_names = FALSE, km = 5)
# for(an in colnames(df)) {
#   decorate_annotation(an, {
#     grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
#   })}

####################
# df <- data.frame(TP53=c(TP53.final),CDKN2A=c(CDKN2A.final), KRAS=c(KRAS.final),Clusters = c(rep("Cluster1", 43), rep("Cluster2",111), rep("Cluster3",30)))
# ha <- HeatmapAnnotation(df = df, col = list(Clusters = c("Cluster1" = "blue4", "Cluster2"="green4", "Cluster3"="red3"), TP53= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white"), CDKN2A= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white"), KRAS= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white")),gap = unit(c(2, 4,4), "mm"))
# Heatmap(Meth.knn.0.3.tumor.NMFcluster,name = "BetaValue", top_annotation = ha, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), cluster_columns = FALSE, cluster_rows = TRUE, show_row_names = FALSE, km = 5, show_column_names = FALSE)
# for(an in colnames(df)) {
#   decorate_annotation(an, {
#     grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
#   })}
################
set.seed(123) ### set seed for getting same cluster in each run
pdf("Heatmap.pdf", width = 10, height = 15)
df <- data.frame(TP53=c(TP53.final),KRAS=c(KRAS.final),CDKN2A=c(CDKN2A.final), SMAD4=c(SMAD4.final),JMY=c(JMY.final),KDM6A=c(KDM6A.final), KMT2C=c(KMT2C.final),KMT2D=c(KMT2D.final),Clusters = c(rep("Cluster1", 43), rep("Cluster2",111), rep("Cluster3",30)))
ha <- HeatmapAnnotation(df = df, col = list(Clusters = c("Cluster1" = "blue4", "Cluster2"="green4", "Cluster3"="red3"), TP53= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white"), 
                                            CDKN2A= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white"), KRAS= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white"), 
                                            KMT2C= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white"), KDM6A= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white"), JMY= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white"),
                                            KMT2D= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white"),SMAD4= c("SNV" = "blue4", "Deletion"="green4", "Amplification"="red3", "NA"="white")), gap = unit(c(5,2,2,2,2,2,2,2), "mm"), show_legend = FALSE)
#show_legend = FALSE if we don't want to print annotation legend
ht <- Heatmap(Meth.knn.0.3.tumor.NMFcluster,name = "BetaValue", top_annotation = ha, col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("blue4", "royalblue", "white", "firebrick2", "darkred")), cluster_columns = FALSE, cluster_rows = TRUE, show_row_names = FALSE, km = 5, show_column_names = FALSE, row_dend_side  = "right", combined_name_fun = NULL)
draw(ht)
for(an in colnames(df)) {
  decorate_annotation(an, {
    grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
  })}
dev.off()
##################
#Get the probes name in each cluster
##################
##Get the probes name in each row clusters 
ht <- Heatmap(Meth.knn.0.3.tumor.NMFcluster,name = "BetaValue", top_annotation = ha, col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("blue4", "royalblue", "white", "firebrick1", "red4")), cluster_columns = FALSE, cluster_rows = TRUE, show_row_names = FALSE, km = 5, show_column_names = FALSE, row_dend_side  = "left", combined_name_fun = NULL)
ht_list <- draw(ht)
row_order(ht_list)
row_order(ht_list)[1] ## This command print the order of row in cluster1
#####################

######################
## Check the SNV, deletion & amplification
tmp <- NMF.cluster.order[RIOK1.sample,]
n <- nrow(tmp)

SNV_Barcode <- unique(substr(Mutect.NonSilent$Tumor_Sample_Barcode, 1, 16))
NMF.cluster_SNV_Barcode <- NMF.cluster.order[SNV_Barcode,]

N <- nrow(NMF.cluster_SNV_Barcode)
m <- length(which(NMF.cluster_SNV_Barcode$cluster==1))
k <- length(which(tmp$cluster==1))
a <- k
b <- n-k
c <- m-k
d <- N+k-n-m
clusterTest <- matrix(c(a,b,c,d), nrow = 2, dimnames = list(c("SNV","NoSNV"),c("InSample", "NotInSample")))
clusterTest
fisher.test(clusterTest, alternative = "g")
fisher.test(clusterTest, alternative = "l")
##################


#Fisher's exact test for deletion
#dim(CNV_all_data_by_genes)
tmp <- NMF.cluster.order[RIOK1.deletion,]
n <- nrow(tmp)
N <- ncol(CNV_all_data_by_genes)
m <- length(which(NMF.cluster.order$cluster==1))
k <- length(which(tmp$cluster==1))
a <- k
b <- n-k
c <- m-k
d <- N+k-n-m
clusterTest <- matrix(c(a,b,c,d), nrow = 2, dimnames = list(c("Deletion","NoDeletion"),c("InSample", "NotInSample")))
clusterTest
fisher.test(clusterTest, alternative = "g")
fisher.test(clusterTest, alternative = "l")

###################
#Fisher's exact test for amplification
tmp <- NMF.cluster.order[SMAD4.amplification,]
n <- nrow(tmp)
N <- ncol(CNV_all_data_by_genes)
m <- length(which(NMF.cluster.order$cluster==1))
k <- length(which(tmp$cluster==1))
a <- k
b <- n-k
c <- m-k
d <- N+k-n-m
clusterTest <- matrix(c(a,b,c,d), nrow = 2, dimnames = list(c("Deletion","NoDeletion"),c("InSample", "NotInSample")))
clusterTest
fisher.test(clusterTest, alternative = "g")
fisher.test(clusterTest, alternative = "l")