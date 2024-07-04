load("BRCA.RData")
BRCA.Clinical.cBioPortal <- read.csv(file = "data_bcr_clinical_data_patient.txt", header = TRUE, sep = "\t", stringsAsFactors=FALSE, check.names = FALSE, row.names=2, skip = 4)
Neetha.freq <- read.csv("BRCA_sample_selection_fusion_frequency.txt", header = TRUE, sep = "\t")
quantile(Neetha.freq$X.Total.Fusions, c(0.2,0.8))
Neetha.freq$quantile <- ifelse(Neetha.freq$X.Total.Fusions <= 2, "Low", ifelse(Neetha.freq$X.Total.Fusions >=154, "Higher", "Middle"))
#Neetha.freq$quantile <- ifelse(Neetha.freq$X.Total.Fusions <= quantile(Neetha.freq$X.Total.Fusions, c(0.2)), "Low", ifelse(Neetha.freq$X.Total.Fusions >=quantile(Neetha.freq$X.Total.Fusions, c(0.8)), "Higher", "Middle"))
## TCGA-A7-A13E1 and TCGA-A7-A26J are two times in "Low"
## Last two times of these two have row name "569" and "640". I have to remove these two from list.
#duplicateID <- c("569","640")
#Neetha.freq1.1 <- Neetha.freq[!rownames(Neetha.freq)%in%duplicateID,]
## Below two are the automatic command to remove duplicated lines. I will take first entry and remove other duplicates from file
dup <- duplicated(Neetha.freq$Sample.ID)
Neetha.freq1.1 <- Neetha.freq[!dup,]
rownames(Neetha.freq1.1) <- Neetha.freq1.1$Sample.ID
Neetha.freq1.cBioPortal <- merge(Neetha.freq1.1, BRCA.Clinical.cBioPortal, by="row.names")
rownames(Neetha.freq1.cBioPortal) <- Neetha.freq1.cBioPortal$Row.names
Neetha.freq1.cBioPortal$Row.names <- NULL
Neetha.freq1.cBioPortal <- Neetha.freq1.cBioPortal[order(Neetha.freq1.cBioPortal$quantile),]


Higher.freq <= Neetha.freq1.cBioPortal[which(Neetha.freq1.cBioPortal$quantile=="Higher"),]
Low.freq <- Neetha.freq1.cBioPortal[which(Neetha.freq1.cBioPortal$quantile=="Low"),]

write.csv(Neetha.freq1.cBioPortal, file = "Neetha.freq1.cBioPortal.csv")

######## Differential expression analysis ######
###################### edgeR #####################
factors <- factor(c(rep("", length(cancerID)), rep("Normal", length(normalID))))
cnts.dgelist <- DGEList(cnts.liver, group=factors)
tmm <- calcNormFactors(cnts.dgelist, method = "TMM")
tmm <- estimateDisp(tmm) ## estimateDisp command do both CommonDisp and TagwiseDis
tmm.DE <- exactTest(tmm)
tmm.DE.top <- topTags(tmm.DE, n=nrow(tmm.DE), sort.by="logFC", p.value = 0.01)
topTags.DEG <- tmm.DE.top[tmm.DE.top$table$FDR < 0.01 & abs(tmm.DE.top$table$logFC) >=1.5,]
diffExp.edgeR <- topTags.DEG$table
diffExp.edgeR$geneID <- sapply(strsplit(rownames(diffExp.edgeR),"\\|"),'[[',2)
diffExp.edgeR$symbol<- sapply(strsplit(rownames(diffExp.edgeR),"\\|"),'[[',1)

save.image(file = "Neetha_BRCA.RData")