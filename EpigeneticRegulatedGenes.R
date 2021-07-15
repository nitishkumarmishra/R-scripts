##### R program for Epigenetic controlled gene
library(pacman)
p_load(ELMER, Homo.sapiens, DMRcate, biomaRt, missMethyl, multtest, impute, matrixStats, ChAMP)
#################################################
getTCGA("CHOL", Meth = TRUE, RNA = TRUE, Clinic = TRUE,   basedir="~", RNAtype = "gene", Methfilter = 0.2)
geneExp <- read.csv("CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", header = TRUE, sep = "\t", row.names = 1)
#geneExp <- geneExp[-1,]
colnames(geneExp) <- gsub("\\.", "-", colnames(geneExp))
save(geneExp, file = "CHOL-firehoseExp.Rda")
load("CHOL-firehoseExp.Rda")
load("CHOL_meth.rda")
colnames(Meth) <- substr(colnames(Meth),1,16)
colnames(geneExp) <- substr(colnames(geneExp),1,16)
Probeinfo <- read.csv("C:/Users/nitish.mishra/Desktop/eMap/TCGA-jhu-usc.edu_TCGA_HumanMethylation450.adf/jhu-usc.edu_TCGA_HumanMethylation450.adf.txt", header = TRUE, sep = "\t", row.names = 1)
#Probeinfo.Meth <- Probeinfo[rownames(Meth),]
#geneSymbol <- intersect(sapply(strsplit(rownames(geneExp),split = "\\|"), '[', 1), Probeinfo.Meth$GENESYMBOL)
###   Make list of common gene symbol
#mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#results <- getBM(attributes = c("entrezgene", "hgnc_symbol"), filters = "hgnc_symbol", values = geneSymbol, mart = mart)
#EntrezID <- sapply(strsplit(rownames(geneExp[geneSymbol,]),split = "\\|"), '[', 2)
##library(mygene)### We can use mygene package for symbol to Ensembel ID, where tt is list of symbol
#queryMany(tt, scopes="symbol", fields="entrezgene", species="human")
############# Preprocess data ###################
#keep <- !(rownames(Meth) %in% ann450k$Name[ann450k$chr %in%c("chrX","chrY")])
#Meth <- Meth[keep,] ## Remove probes from chrX and chrY
Meth <- rmSNPandCH(Meth, dist = 10, mafcut = 0.05, rmXY = TRUE, rmcrosshyb = TRUE)
Meth <- Meth[-grep("rs", rownames(Meth)),]## Remove rs probes
Tumor.Meth <- Meth[,which(substr(colnames(Meth),14,15)=="01")]
Normal.Meth <- Meth[,which(substr(colnames(Meth),14,15)=="11")]
numNAs <- rowSums(is.na(Tumor.Meth))
Tumor.Meth <- Tumor.Meth[!(numNAs > dim(Tumor.Meth)[2]*0.25),]## Remove probe which have >25% NAs
numNAs <- rowSums(is.na(Normal.Meth))
Normal.Meth <- Normal.Meth[!(numNAs > dim(Normal.Meth)[2]*0.25),]
commonID <- intersect(rownames(Tumor.Meth), rownames(Normal.Meth))
Tumor.Meth <- Tumor.Meth[commonID,]
Normal.Meth <- Normal.Meth[commonID,]
Tumor.impue <- impute.knn(Tumor.Meth, k = 15, rowmax = 0.25, colmax = 0.25, rng.seed = 12345)
Tumor.impue <- Tumor.impue$data
Normal.impue <- impute.knn(Normal.Meth, k = 15, rowmax = 0.25, colmax = 0.25, rng.seed = 12345)
Normal.impue <- Normal.impue$data
myNorm <- champ.norm(beta = Tumor.impue, fromIDAT = F, methValue = "B", norm = "BMIQ", fromFile = FALSE, betaFile = Tumor.impue, filterXY = TRUE, plotBMIQ = TRUE, resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
Tumor.impue.BMIQ <- myNorm$beta
myNorm <- champ.norm(beta = Normal.impue, fromIDAT = F, methValue = "B", norm = "BMIQ", fromFile = FALSE, betaFile = Normal.impue, filterXY = TRUE, plotBMIQ = TRUE, resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
Normal.impue.BMIQ <- myNorm$beta
Tumor.exp <- geneExp[,colnames(Tumor.Meth)]
Normal.exp <- geneExp[,colnames(Normal.Meth)]
Probeinfo.merge <- merge(Probeinfo, probe.features, by = "row.names") ## Merge chAMP and TCGA DNA annotation file
###### Remove probes which have mean beta value >0.3 in normal 
keep <- !rowMeans(Normal.impue.BMIQ) >0.3
Tumor.impue.BMIQ.keep <- Tumor.impue.BMIQ[keep,]
Tumor.hyper <- (Tumor.impue.BMIQ.keep >=0.3)+0 ##convert matrix in binary at beta >= 0.3
numNAs <- rowSums(Tumor.hyper)
### probe at least 5% hypermethylated
Hyper <- Tumor.impue.BMIQ.keep[(numNAs > round(dim(Tumor.impue.BMIQ.keep)[2]*0.05)),]
#Hyper.bin <- (Hyper >=0.3)+0
################ Function to call ENTREZ-GENE-ID and get 1.5kb up and downstream probes from TSS ########
symbol <- sapply(strsplit(rownames(geneExp),split = "\\|"), '[', 1)
id <- sapply(strsplit(rownames(geneExp),split = "\\|"), '[', 2)
#aux <- strsplit(row.names(exp),"\\|")
#Gene_Symbol  <- unlist(lapply(aux,function(x) x[1]))
tmp <- cbind(symbol, id)
symbol.entrezID <- as.data.frame(tmp[-c(1:29),])### Remove first 29 line which have no HGNC symbol
rownames(symbol.entrezID) <- symbol.entrezID$id
#####################################################
#use Entrez GeneID as a string, not numeric or factor
library(Homo.sapiens)
keytypes(Homo.sapiens)
txs <- transcriptsBy(Homo.sapiens, 'gene', col='GENEID')
library(FDb.InfiniumMethylation.hg19)
hm450 <- get450k()
#use Entrez GeneID as a string, not numeric or factor
getProbes<-function(geneID){
  temp<-txs[[geneID]]
  if(is.null(temp)){
    probes<-NULL
  }else{
    ### Use 250bp from TSS suggested by Duttke et.al. Mol. Cell (2015) 57, 674-684
    upstream.probes<-names(subsetByOverlaps(hm450,flank(temp,250)))
    downstream.probes<-names(subsetByOverlaps(hm450,flank(temp,-250,start=TRUE)))
    probes<-unique(c(upstream.probes,downstream.probes))
  }
  return(probes)
}
#example
#getProbes('1234')
############################ Final line to get list of probes for all genes ##########
y <- as.character(symbol.entrezID$id)
Probe250kb <- lapply(y, getProbes)
save.image("EpigeneticRegulatedGenes.RData")
#tmp <- do.call(rbind, Probe250kb)## It will convert in matrix of same row size. Any row have less element then it will filled by same element repeatedly
tt <- which(Probe250kb=="NULL")## List with of index in "y" with no TSS information 
tt1 <- which(lengths(Probe250kb)==0)##Total with NA and character(0)
tt2 <- setdiff(tt1, tt)## This is the index of genes in "y" which have no probes in 250bp region
TSS250 <- as.data.frame(do.call(rbind, Probe250kb))## It will convert in matrix of same row size. Any row have less element then it will filled by same element repeatedly
## In TSS250 NULL and character(0) rows are absent
rownames(TSS250)<-y[(which(lengths(Probe250kb)!=0))] ## EntrezID in rowname of TSS250
TSS250$id <- rownames(TSS250)
TSS250 <- merge(TSS250, symbol.entrezID, by="id")
TSS250.ID <- paste(as.character(TSS250$symbol), TSS250$id, sep = "|")
rownames(TSS250) <- TSS250.ID 
TSS250.1 <- TSS250; TSS250.1$id <- NULL; TSS250.1$symbol <- NULL
matrix <- as.matrix(TSS250.1)
matrix1 <- t(matrix)
list <- apply(matrix1,2, unique)
d1 <- data.frame(geneID = rep(names(list), sapply(list, length)), Probe = unlist(list))## This will print gene name and probe in each row, if gene have more than one probe then it will print more than one line
##############################################
geneExp.Hyper <- geneExp[rownames(matrix),colnames(Hyper)]
geneExp.Hyper.log <- log2(geneExp.Hyper+1)
###############################################
d2 <- as.data.frame(Hyper)
d2$Probe <- rownames(d2)
d2 <- merge(d1, d2, by="Probe")
##############################################
t1 <- geneExp.Hyper.log
t1$geneID <- rownames(t1)
d3 <- merge(d2, t1, by="geneID")
#apply(Hyper, 1, function(x) mean(x[x< 0.3]))## Mean of rows with Hyper beta less than 0.3
#apply(Hyper,1, function(x) which(x < 0.3)) ## This command will make Hyper.bin type list
meth <- names(d3[3:38])
exp <- names(d3[39:74])
identical(substring(meth, 1, 15), substring(exp, 1, 15))
meth <- as.matrix(d3[,3:38])
exp <- as.matrix(d3[,39:74])

hypo <- apply(meth,1, function(x) which(x < 0.3))## Substitute of Hyper.bin
hyper <- apply(meth,1, function(x) which(x >= 0.3))

a1 <- matrix(nrow=14157, ncol=36)
a2 <- matrix(nrow=14157, ncol=36)
for(i in 1:14157)
{
  a <- unlist(hypo[i])
  a1[i,a] <- exp[i,a]
  b <- unlist(hyper[i])
  a2[i,b] <- exp[i,b]
}
d3$hyperExpMean <- rowMeans(a2, na.rm=TRUE)
d3$hypoExpMean <- rowMeans(a1, na.rm=TRUE)
d3$hypoExpSD <- rowSds(a1, na.rm=TRUE)
d3$hyperExpSD <- rowSds(a2, na.rm=TRUE)
d3$quartile10.hypo <- apply(a1, 1, quantile, probs=0.1, na.rm=TRUE)
d3$lower10SD1.282 <- d3$hypoExpMean+(-1.282*d3$hyperExpSD)
d3$totalHyper <- apply(a2, 1, function(x) length(which(!is.na(x))))
ss <- subset(a2 < d3$hypoExpMean) ## Probes corresponding gene in Hyper have expression be;ow mean of Hypo group
d3$totalHyperWithGreaterMeanHypo <- rowCounts(ss, na.rm = TRUE, value = TRUE)
ss <- subset(a2 > d3$hypoExpMean) ## Probes corresponding gene in Hyper have expression be;ow mean of Hypo group
d3$totalHyperWithLowerMeanHypo <- rowCounts(ss, na.rm = TRUE, value = TRUE)
d3$percentageHyper <- (d3$totalHyperWithGreaterMeanHypo/d3$totalHyper)*100
########################################################################
a3 <- matrix(nrow=14157, ncol=36)
a4 <- matrix(nrow=14157, ncol=36)
for(i in 1:14157)
{
  a <- unlist(hypo[i])
  a3[i,a] <- meth[i,a]
  b <- unlist(hyper[i])
  a4[i,b] <- meth[i,b]
}
###### calculate correlation between Hyper methylated probe beta value and corresponding gene's expression
ll <- mapply(function(x,y)cor.test(a2[x,],a4[y,], method = "pearson", alternative = "l"),
             1:nrow(a2),
             1:nrow(a4),
             SIMPLIFY=FALSE)
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
adjust <- mt.rawp2adjp(p.value, proc=(c("Bonferroni")))
adjust.order <- adjust$adjp[order(adjust$index),]
##################################
d3$corr <- cor.value
d4 <- d3[(which(d3$hyperExpMean < d3$lower10SD1.282)),]
d5 <- d3[(which(d3$hyperExpMean < d3$lower10SD1.282  & d3$corr <= -0.25)),]
############################################################################
#qq <- table(d5$geneID)
qq <- table(d5$geneID)
qq1 <- table(d3$geneID)
qq2 <- cbind(qq, qq1)
qq3 <- qq2[which(qq2[,1] >=1),]
EpigeneticSlencedGene <- qq3[which(qq3[,1] > qq3[,2]/2),]
detail.EpigeneticSlencedGene <- d5[d5$geneID%in%rownames(EpigeneticSlencedGene),]
############################################################################
EpigeneticGene <- unique(detail.EpigeneticSlencedGene$geneID)
sapply(strsplit(as.character(EpigeneticGene), split = "\\|"), '[', 1)