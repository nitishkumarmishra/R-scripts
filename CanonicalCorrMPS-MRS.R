## This program calculate canonical correlation between expression and mehylation
## This program developed by using getNearestTSS command from FDb.InfiniumMethylation.hg19
## Unlike Hui Shen program it's working for all gene, even if there is no CpG in genomica range of given gene
library('org.Hs.eg.db')
library(Homo.sapiens)
library(biomaRt)
library(FDb.InfiniumMethylation.hg19)
library(CCA)
#######################################################
Expr <- readRDS("Expr.me.rds")
genes <- rownames(Expr)
Meth <- readRDS("Meth.mm.rds")
#genes <- read.csv("TSS-HGNC-name.csv") ### TSS-HGNC.CSV is the list of TCGA genes in RNASeqV2 level-3 data
#yy<-as.character(genes$hgnc, na.rm=FALSE)
keytypes(Homo.sapiens)
hm450 <- get450k()
TSS.nearest <- getNearestTSS(hm450)##  Make file for nearest TSS for each CpGs
### Meth have selected CpG's, we have to select probes which are in data and within 1.5Kb from TSS
### Keep common CpG's only
CpGlist <- rownames(Meth)
CpGlist.TSS1500 <- rownames(TSS.nearest[TSS.nearest$distance <=1500,])
CpGlist <- intersect(CpGlist, CpGlist.TSS1500)
TSS.nearest <- TSS.nearest[CpGlist,]
TSS.nearest.uniq.gene <- unique(TSS.nearest$nearestGeneSymbol)
gene <- intersect(genes, TSS.nearest.uniq.gene)
#gene <- intersect(yy, TSS.nearest.uniq.gene)
#tmp <- TSS.nearest[grep("^A1CF$", TSS.nearest$nearestGeneSymbol),]
#rownames(tmp[tmp$distance <=1500,])
getProbes<-function(geneID){
  #tmp <- TSS.nearest[grep(geneID, TSS.nearest$nearestGeneSymbol, fixed = TRUE),]## This command have problem in NAT1, NAT10 and NAT14
  # In case of NAT1 it give list of all NAT1 (NAT1, NAT10, NAT14 etc)
  tmp <- TSS.nearest[which(TSS.nearest$nearestGeneSymbol==geneID),]
  if(is.null(tmp)){
    probes <- NULL}
  else{
    probes <- rownames(tmp[tmp$distance <=1500,])
  }
  return(probes)
}
#listprobes <- sapply(gene1, getProbes)
#listprobes <- sapply(gene, getProbes)
########################################################
#### Code for extracing beta value and expression value
####### Get expression value for given gene name ######
getExpr<-function(geneID){
  tmp <- Expr[which(rownames(Expr)==geneID),]
  return(tmp)
}
####### Get beta value for given CpG name ######
getMeth<-function(probeID){
tmp <- Meth[which(rownames(Meth)==probeID),]
return(tmp)
}
####### Get methylation for each CpG in promoter of given gene gene (1.5 Kb +/- from TSS) ######
getMethVal<-function(geneID){
  probes <- unlist(sapply(geneID, getProbes))
  probes <- intersect(probes, CpGlist)
  meth <- sapply(probes, getMeth)
  return(meth)
}
#can <- cancor(t(matrix(unlist(tmp1), ncol = ncol(Expr))),t(matrix(unlist(tmp), ncol = ncol(Meth))))
####### Get canonical correlation ######
getCanCorr <- function(geneID){
  meth <- getMethVal(geneID)
  #meth <- t(matrix(unlist(meth), ncol = ncol(Meth)))
  #expr <- sapply(geneID, getExpr)
  expr <- getExpr(geneID)
  expr <- t(t(expr))
  colnames(expr) <- geneID
  #expr <- t(matrix(unlist(expr), ncol = ncol(Expr)))
  can <- cc(meth, expr)
  return(can)
}
#lapply(gene2, getCanCorr)### Here lapply make list but sapply is making matrix which is wrong

####### Get $xcoef for each CpG ######
getXcoef <- function(geneID){
  can <- getCanCorr(geneID)
  xcoef <- can$xcoef
  return(xcoef)
}
####### Get MPS (Methylation pattern score) for each gene ####
## We calculate MPS for each gene in each sample####
getMPS <- function(geneID){
  meth.gene <- getMethVal(geneID)
  xcoef <- getXcoef(geneID)
  mps <- meth.gene%*%xcoef
  colnames(mps) <- geneID
  return(mps)
}
### Test for fisr 25 genes ###########
gene2 <- gene[1:25]
MPS <- t(sapply(gene2, getMPS))### This is the MPS, but colname is missing
colnames(MPS) <- colnames(Meth) ## Manually add colname in mtrix
med <- cbind(rowMeans(MPS),rowMedians(MPS))## It didn't give rowname (gene name)
#apply(MPS, 1, median) # This command give rowname
colnames(med) <- c("Mean","Median")

