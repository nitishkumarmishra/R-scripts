## BRCA_Samples_fusion_list.txt is tab separated (three column) file. I edit file in notepad++ before using in R.
## Second column is TCGA id and third column is comma (,) separated fusion genes list
setwd("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/GDC-Harmonized/BRCA")
load("BRCA_Neetha_GSEA_Analysis.RData")
rm(list= ls()[!(ls() %in% c('BRCA.PAM50', 'Neetha.freq','BRCA.Clinical','BRCA.Clinical.cBioPortal', 'BRCA.FPKM', 'BRCA.FPKM.UQ', 'BRCA.HTSeq.Counts', 'TCGA.Harmonized.RNASeq.Annotation'))])
df <- read.csv("BRCA_Samples_fusion_list.txt", header = FALSE, sep = '\t', stringsAsFactors = FALSE)
s <- strsplit(df$V3, split = ",")
s <- data.frame(V1 = rep(df$V2, sapply(s, length)), V2 = unlist(s), stringsAsFactors = FALSE)

# tx <- table(s$V2); tx <- tx[order(tx, decreasing = TRUE)]
# tx1 <- table(s$V1); tx1 <- tx[order(tx1, decreasing = TRUE)]
#uniq_genes <- unique(unlist(strsplit(s$V2, split = "&", fixed = TRUE)))
# t1 <- aggregate(V2~., s, FUN=toString); t1.1 <- aggregate(V1~., s, FUN=toString)
# library(stringr)
# count <- str_count(t1.1$V1, ', ') ## count number of comma (,) i.e. number of sample for fusion genes

s1 <- strsplit(s$V2, split = "&", fixed = TRUE) ## or strsplit(s$V2,"\\&")# both are same
s1.1 <- sapply(strsplit(s$V2,"\\&"),'[[',1)
s1.2 <- sapply(strsplit(s$V2,"\\&"),'[[',2)
s1 <- data.frame(TCGAID = s$V1, Gene1 = s1.1, Gene2 = s1.2)

t2.1 <- aggregate(TCGAID~Gene1, s1, FUN=toString)
colnames(t2.1) <- gsub("Gene1", "Gene", colnames(t2.1))
t2.2 <- aggregate(TCGAID~Gene2, s1, FUN=toString)
colnames(t2.2) <- gsub("Gene2", "Gene", colnames(t2.2))
t2.1$TCGAID <- gsub(", ", ",", t2.1$TCGAID); t2.2$TCGAID <- gsub(", ", ",", t2.2$TCGAID)
pp <- merge(t2.1, t2.2, by = "Gene")
colnames(pp) <- gsub("Gene.x", "Gene", colnames(pp))
TCGAID <- paste(pp$TCGAID.x,pp$TCGAID.y,sep=",")
pp <- cbind(pp, TCGAID)
#pp <- subset(pp, select = -c(Row.names, TCGAID.x, Gene.y, TCGAID.y )) # or  pp <- pp[,c("Gene", "TCGAID")]
pp <- subset(pp, select = -c(TCGAID.x, TCGAID.y )) # or  pp <- pp[,c("Gene", "TCGAID")]

pp.1 <- t2.1[!t2.1$Gene %in% pp$Gene,]
pp.2 <- t2.2[!t2.2$Gene %in% pp$Gene,]
rownames(pp.1) <- pp.1$Gene; rownames(pp.2) <- pp.2$Gene; rownames(pp) <- pp$Gene

m <- rbind(pp.1, pp.2, pp); m <-  m[order(rownames(m)),]
m.1 <- strsplit(m$TCGAID, ","); names(m.1) <- m$Gene
m.1 <- lapply(m.1, unique)
len <- sapply(m.1, length); nm <- names(len[len>=3]) ## Name of genes (fusion gene) which are present in 3 or more samples 
m.1.nm <- m.1[names(m.1)%in% nm]

TCGA.Harmonized.RNASeq.Annotation <- read.csv("TCGA_GDC_Harmonided_RNASeq-GFT.txt", header = TRUE, sep = "\t")
rownames(TCGA.Harmonized.RNASeq.Annotation) <- TCGA.Harmonized.RNASeq.Annotation$ENSB_ID
TCGA.Harmonized.RNASeq.Annotation.subset <- subset(TCGA.Harmonized.RNASeq.Annotation, select=-c(Chr, Reference, Start,End, Strand, ENSB_ID)) 
proteinCoding <- merge(log2(BRCA.HTSeq.Counts+1), TCGA.Harmonized.RNASeq.Annotation.subset, by="row.names")
exp.proteinCoding <- proteinCoding[which(proteinCoding$Gene_status=="protein_coding"),]
rownames(exp.proteinCoding) <- exp.proteinCoding$Row.names 
exp.proteinCoding$Row.names <- NULL; exp.proteinCoding$Gene_status <- NULL
exp.proteinCoding.nm <- exp.proteinCoding[exp.proteinCoding$Gene_Symbol%in%nm,]
uniq.TCGA <- unique(unlist(m.1.nm))
#uniq.TCGA <- unique(gsub(" ", "", uniq.TCGA)) # Some TCGA id have space (" ") before name, I have to remove it.
exp.proteinCoding.nm.TCGA <- exp.proteinCoding.nm[,substr(colnames(exp.proteinCoding.nm), 1, 12)%in%uniq.TCGA] 
exp.proteinCoding.nm.TCGA <- merge(exp.proteinCoding.nm.TCGA, TCGA.Harmonized.RNASeq.Annotation.subset, by = "row.names")
# gene 'RPL41', 'SERPINA3', 'SLC25A6' have two ENSG id. But one ENSG of SLC25A6 have all zero's (line number 246). So I will remove line number 246.
# Similarly 245 ENSG00000279483.1 (245) and ENSG00000273259.2 (242) are also uncharacterized. So I removed these two from data.
exp.proteinCoding.nm.TCGA <- exp.proteinCoding.nm.TCGA[-c(242,245,246),]
rownames(exp.proteinCoding.nm.TCGA) <- exp.proteinCoding.nm.TCGA$Gene_Symbol
exp.proteinCoding.nm.TCGA <- subset(exp.proteinCoding.nm.TCGA, select = -c(Row.names, Gene_status, Gene_Symbol))
#####################################################################
mean.fusion <- vector("numeric", length = nrow(exp.proteinCoding.nm.TCGA))
mean.nonfusion <- vector("numeric", length = nrow(exp.proteinCoding.nm.TCGA))
p.value <- vector("numeric", length = nrow(exp.proteinCoding.nm.TCGA))
t.statistic<- vector("numeric", length = nrow(exp.proteinCoding.nm.TCGA))
p.value.var <- vector("numeric", length = nrow(exp.proteinCoding.nm.TCGA))
F.statistic.var<- vector("numeric", length = nrow(exp.proteinCoding.nm.TCGA))

for(i in 1:nrow(exp.proteinCoding.nm.TCGA)){
  fusion <- which(substr(colnames(exp.proteinCoding.nm.TCGA), 1, 12)%in%m.1.nm[[i]])
  nonfusion <- which(!substr(colnames(exp.proteinCoding.nm.TCGA), 1, 12)%in%m.1.nm[[i]])
  tmp <- t.test(exp.proteinCoding.nm.TCGA[i,fusion], exp.proteinCoding.nm.TCGA[i, nonfusion])
  tmp1 <- var.test(as.numeric(exp.proteinCoding.nm.TCGA[i,fusion]), as.numeric(exp.proteinCoding.nm.TCGA[i, nonfusion]))
  p.value[i] <- tmp$p.value
  t.statistic[i] <- tmp$statistic
  mean.fusion[i] <- mean(as.numeric(exp.proteinCoding.nm.TCGA[i,fusion]))
  mean.nonfusion[i] <- mean(as.numeric(exp.proteinCoding.nm.TCGA[i, nonfusion]))
  p.value.var[i] <- tmp1$p.value
  F.statistic.var[i] <- tmp1$statistic
}
Fusion.test <- as.data.frame(cbind(mean.fusion, mean.nonfusion, t.statistic, p.value, F.statistic.var,p.value.var))
Fusion.test$p.value.fdr <- p.adjust(Fusion.test$p.value, method = "fdr")
Fusion.test$p.value.var.fdr <- p.adjust(Fusion.test$p.value.var, method = "fdr")
Fusion.test <- subset(Fusion.test, select = c(mean.fusion, mean.nonfusion, t.statistic, p.value, p.value.fdr, F.statistic.var,p.value.var, p.value.var.fdr))
rownames(Fusion.test) <- rownames(exp.proteinCoding.nm.TCGA)

################## Remove unnecessay extra files ####################
rm(list= ls()[!(ls() %in% c('Fusion.test','exp.proteinCoding.nm.TCGA','exp.proteinCoding','m.1','m.1.nm','nm','BRCA.PAM50', 'Neetha.freq','BRCA.Clinical','BRCA.Clinical.cBioPortal', 'BRCA.FPKM', 'BRCA.FPKM.UQ', 'BRCA.HTSeq.Counts', 'TCGA.Harmonized.RNASeq.Annotation'))])
#####################################################################
Fusion.test <- na.omit(Fusion.test)
Fusion.test[Fusion.test$p.value.fdr < 0.05,]
################### save data #######################################
save.image(file = "Comparative_Fusion.RData")
#####################################################################
### upper quantile normalization 
data.mat <- data.mat[rowSums(data.mat) > 0,] ## Only retain rows that have a positive total expression per gene
data.quantileAll <- apply(data.mat, 2, function(x){quantile(x, 0.75)})
data.quantileExpressed <-   apply(data.mat, 2, function(x){quantile(x[x>0], 0.75)})
data.norm <- t(t(data.mat) / data.quantileExpressed)
norm_edata = preprocessCore::normalize.quantiles(as.matrix(data.mat)) ## Quantile normalization by using preProcessCore
