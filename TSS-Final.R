library(matrixStats)
library(samr)
library(wateRmelon)
library(data.table) #----When CSV file is big then we can use fread rather than read.csv -----#

x <- as.matrix(read.csv(file="MethylationBetaTable.csv", header=TRUE, sep = ",", row.names = 1, as.is=TRUE))

#----- BELOW STEP ONLY WHEN CSV FILE IS BIG AND CREATING PROBLEM WITH read.csv/read.table -----------#
x <- fread('MethylationBetaTable.txt', sep="auto", sep2="auto", nrows=-1, header = TRUE, drop=196)
a <- x$V1 #####(This is the list of rownames)
x$V1 <-NULL
rownames(x) <- a ### Now this x is data.frame 
matrix <- data.matrix(x) ### Finally this is the "t" matrix
#rm(list = c('a','x'))
#---------------------------------------------------------------------------------------------------#

HM450ProbeInfo <- as.matrix(read.table(file="HM450ProbeInfo.txt", header=TRUE, row.names = 1, as.is=TRUE)) # HM450ProbeInfo.txt is file of HM450 annotation
#HM450ProbeInfo.txt is subset of jhu-usc.edu_TCGA_HumanMethylation450/IMA which have only 3 column "UCSC_REFGENE_GROUP", "UCSC_CPG_ISLANDS_NAME", RELATION_TO_UCSC_CPG_ISLAND"
TSS200 <- rownames(HM450ProbeInfo[grep("TSS200", HM450ProbeInfo), ])## This command will print probes which have TSS200
TSS1500 <- rownames(HM450ProbeInfo[grep("TSS1500", HM450ProbeInfo), ])
TSSaLL <- rownames(HM450ProbeInfo[grep("TSS", HM450ProbeInfo), ])
UTR5 <- rownames(HM450ProbeInfo[grep("5'UTR", HM450ProbeInfo), ])## This command will print probes which have TSS200
UTR3 <- rownames(HM450ProbeInfo[grep("3'UTR", HM450ProbeInfo), ])
Body <- rownames(HM450ProbeInfo[grep("Body", HM450ProbeInfo), ])
Exon1st <- rownames(HM450ProbeInfo[grep("1stExon", HM450ProbeInfo), ])
ProbesAllSixClass <- unique(c(TSS200, TSS1500, UTR3, UTR5, Body, Exon1st))## Unique list of probes from all six class (TSS*, UTR*, Exon1st, Body)
#TSS1.5KbUPdOwn <- unlist(sapply(y, getProbes))## List of probes in 1.5 Kb Up and downstream of TSS using GenomicFeatures

#######################################################################
t <- matrix[rownames(matrix) %in% ProbesAllSixClass,]
a <- subset(t, rowSds(t[,1:184]) <= 0.15) ## select probe which have SD less or equal to 0.15 i.e total 48887 probes
b <- subset(t, rowSds(t[,185:194]) <= 0.15)## select probe which have SD <= 0.15 only in tumor samples
m <- row.names(a)
n <- row.names(b)
p <- intersect(m,n)
afterSd <- t[rownames(t) %in% p,]
a1 <- afterSd[,1:184]
b1 <- afterSd[,185:194]
i <- (cbind(rowMeans(a1), rowSds(a1), rowMedians(a1)))
j <- (cbind(rowMeans(b1), rowSds(b1), rowMedians(b1)))
z <- cbind(i,j)
z1 <- z[abs((z[,3] - z[,6])) >= 0.1,]#dim(z1) = 846
m1 <- row.names(z1)
t1 <- afterSd[rownames(afterSd) %in% m1,]
###t1[rownames(t1)==p[],]
##################################################
###beta2m is lumi function called in wateRMelon
#system.time(beta.norm <- betaqn(t1))
M.norm <- beta2m(t1)
y <- c(rep(2,184),rep(1,10))
data <- list(x=t1,y=y,logged2=FALSE,genenames=paste(row.names(t1)), geneid=paste(row.names(t1)))
#data <- list(x=M.norm,y=y,logged2=TRUE,genenames=paste(row.names(M.norm)), geneid=paste(row.names(M.norm)))

###################
samr.obj <- samr(data, resp.type = "Two class unpaired", nperms = 100, random.seed = 123456)
delta.table <- samr.compute.delta.table(samr.obj, nvals=50)
del <- 0.1
siggenes.table <- samr.compute.siggenes.table(samr.obj, del, data, delta.table, min.foldchange = 1.5, all.genes=FALSE)
###############
up <- siggenes.table$genes.up
down <- siggenes.table$genes.lo
upDown <- rbind(up,down)
low <- upDown[as.numeric(upDown[,8])<1,] ### Numerator in samr table is difference me mean of two group
head(low[order(low[,7], decreasing=FALSE),])## sort on the basis of column 7
