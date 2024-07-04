library(edgeR)
library(limma)
library(DESeq2)
setwd("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/GDC-Harmonized/BRCA")
load("BRCA.RData") ### Load TCGA harmonized data for BRCA from oneDrive
####### Prepairing file for fusion frequency and gene expression analysis #####
BRCA.Clinical.cBioPortal <- read.csv(file = "data_bcr_clinical_data_patient.txt", header = TRUE, sep = "\t", stringsAsFactors=FALSE, check.names = FALSE, row.names=2, skip = 4)
TCGA.Harmonized.RNASeq.Annotation <- read.csv("TCGA_GDC_Harmonided_RNASeq-GFT.txt", header = TRUE, sep = "\t")
rownames(TCGA.Harmonized.RNASeq.Annotation) <- TCGA.Harmonized.RNASeq.Annotation$ENSB_ID
TCGA.Harmonized.RNASeq.Annotation.subset <- subset(TCGA.Harmonized.RNASeq.Annotation, select=-c(Chr, Reference, Start,End, Strand, ENSB_ID)) 
BRCA.PAM50 <- read.csv(file = "PAM50_classification_cell_Ciriello etal_2015.txt", header = TRUE, sep = "\t", stringsAsFactors=FALSE, check.names = FALSE, row.names=1, skip = 2)
Neetha.freq <- read.csv(file = "Updated_recurrent_fusion_sample_list.txt", header = TRUE, sep = "\t", stringsAsFactors=FALSE, check.names = FALSE)
dup <- duplicated(Neetha.freq$`TCGA ID`)
Neetha.freq <- Neetha.freq[!dup,]
rownames(Neetha.freq) <- Neetha.freq$`TCGA ID`

Neetha.freq.1 <- Neetha.freq[which(Neetha.freq$`Number of fusions` < 250),]
Neetha.freq.1$quantile <- ifelse(Neetha.freq.1$`Number of fusions` <= quantile(Neetha.freq.1$`Number of fusions`, c(0.25)), "Low", ifelse(Neetha.freq.1$`Number of fusions` >=quantile(Neetha.freq.1$`Number of fusions`, c(0.75)), "Higher", "Middle"))
Neetha.freq.1.cBioPortal <- merge(Neetha.freq.1[,c("Number of fusions", "quantile")], BRCA.Clinical.cBioPortal, by="row.names")

#### write CSV file ############
write.csv(Neetha.freq.1.cBioPortal, file = "neetaCbioPortal.csv")
##### Make expression file for DEG analysis #############
exp.High <- BRCA.HTSeq.Counts[,substr(colnames(BRCA.HTSeq.Counts), 1, 12)%in%rownames(Neetha.freq.1[grep("Higher", Neetha.freq.1$quantile),])]
exp.Low <- BRCA.HTSeq.Counts[,substr(colnames(BRCA.HTSeq.Counts), 1, 12)%in%rownames(Neetha.freq.1[grep("Low", Neetha.freq.1$quantile),])]
exp.High <- exp.High[,substr(colnames(exp.High), 13, 15)=="-01"]
exp.Low <- exp.Low[,substr(colnames(exp.Low), 13, 15)=="-01"]
exp.fusion <- cbind(exp.High, exp.Low)

######## Filtering of genes #####
exp.fusion <- exp.fusion[rowSums(exp.fusion==0)< ncol(exp.fusion)*0.2,] ## Remove all gene which have 25% zero's
keep <- rowSums(cpm(exp.fusion)>1) >= ncol(exp.fusion)*0.5 #### Select only genes which have have CPM > 1 for >=50% samples
exp.fusion <- exp.fusion[keep,]

exp.fusion.proteinCoding <- merge(exp.fusion, TCGA.Harmonized.RNASeq.Annotation.subset, by="row.names")
exp.fusion.proteinCoding <- exp.fusion.proteinCoding[which(exp.fusion.proteinCoding$Gene_status=="protein_coding"),]
rownames(exp.fusion.proteinCoding) <- exp.fusion.proteinCoding$Row.names ## Several gene have more than 1 protein coding transcrip
# Here I still use EnsembleId later on after DEG analysis I will use HGNC gene symbol
exp.fusion.proteinCoding <- subset(exp.fusion.proteinCoding, select=-c(Row.names, Gene_status, Gene_Symbol))

######## Differential expression analysis ######
###################### edgeR #####################
factors <- factor(c(rep("High", length(exp.High)), rep("Low", length(exp.Low))))
cnts.dgelist <- DGEList(exp.fusion.proteinCoding, group=factors)
tmm <- calcNormFactors(cnts.dgelist, method = "TMM")
design <- model.matrix(~0 + factors)
colnames(design) <- c("High", "Low")
tmm <- estimateDisp(tmm, design = design) ## estimateDisp command do both CommonDisp and TagwiseDis
tmm.DE <- exactTest(tmm, pair = c("Low", "High"), prior.count = 0.1) # by default second factor vs first 
tmm.DE <- topTags(tmm.DE, n=nrow(tmm.DE), sort.by="logFC", adjust.method = "BH")
topTags.DEG <- tmm.DE[tmm.DE$table$FDR < 0.05 & abs(tmm.DE$table$logFC) >=1,]
diffExp.edgeR <- topTags.DEG$table

merge(diffExp.edgeR, TCGA.Harmonized.RNASeq.Annotation.subset, by="row.names")
#################### DESeq2 ##########################
cond <- factor(c(rep(2, length(exp.High)), rep(1, length(exp.Low))))
dds <- DESeqDataSetFromMatrix(exp.fusion.proteinCoding, DataFrame(cond), ~ cond)
dds <- DESeq(dds)
resmiRNA <- results(dds, pAdjustMethod = "BH", contrast = c("cond", "2", "1"), lfcThreshold = 0, alpha = 0.05, altHypothesis = "greaterAbs") ## Here I have to define coefficients
## cond is design, 2 is high and 1 is low. This is analysis for 2 Vs 1 (High Vs Low).
resmiRNAOrdered <- resmiRNA[order(resmiRNA$log2FoldChange, decreasing = TRUE),]
resmiRNASig <- subset(resmiRNAOrdered, (padj < 0.05& abs(log2FoldChange) >=1))
diffExp.DESeq2 <- as.data.frame(resmiRNASig)

list_of_data = list(diffExp.edgeR , diffExp.DESeq2)
common_names = Reduce(intersect, lapply(list_of_data, row.names))

rm(dup, cont.matrix, design, factors, adj.method, adj.pval, logFC, raw.pval, cond)

save.image(file = "BRCA_Neetha_GSEA_Analysis.RData")

## This part is extra, I just keep it.
## I didn't find DEG in Limma+Voom
############### Limma + voom #####################
design <- model.matrix(~0 + factor(c(rep("High", length(exp.High)), rep("Low", length(exp.Low)))))
colnames(design) <- c("High", "Low")
# Below two line are same. Here highr value (1 in case of 0,1 and 2 in case of 1,2) is second line of design.
#design <- model.matrix(~0 + factor(c(rep(1, length(exp.High)), rep(0, length(exp.Low)))))
#colnames(design) <- c("Low", "High")
cont.matrix <- makeContrasts("High-Low", levels = design)
factors <- factor(c(rep("High", length(exp.High)), rep("Low", length(exp.Low))))
cnts.dgelist <- DGEList(exp.fusion.proteinCoding, group=factors)
tmm <- calcNormFactors(cnts.dgelist, method = "TMM")
v <- voom(tmm, design, plot = TRUE)
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
adj.method = "BH"; adj.pval = 0.05; raw.pval = 0.01; logFC = 1
aradeger <- topTable(fit2, adjust.method = adj.method, number = nrow(exp.fusion), sort.by='logFC', coef = "High-Low")
diffExp.voom <- data.frame(aradeger[aradeger$adj.P.Val < adj.pval & abs(aradeger$logFC) >=logFC,])


