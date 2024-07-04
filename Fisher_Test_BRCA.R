## Fisher's exact test for HER2+,ER+ and other clinical features
# We used recurrent fusion data. For this study we follow two strategy==
#1::By using recurrent fusion data and whole TCGA clinical data as suggested in IPA and other pathway tools
#2::By using only recurrent fusion data, as suggested in Wikipedia and other sources.
setwd("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/GDC-Harmonized/BRCA")
load("Comparative_Fusion.RData")

BRCA.Clinical.cBioPortal.select <- subset(BRCA.Clinical.cBioPortal, select = c(ER_STATUS_BY_IHC, PR_STATUS_BY_IHC, IHC_HER2))
# df <- read.csv("BRCA_Samples_fusion_list.txt", header = FALSE, sep = '\t', stringsAsFactors = FALSE)
# s <- strsplit(df$V3, split = ",")
# s <- data.frame(V1 = rep(df$V2, sapply(s, length)), V2 = unlist(s), stringsAsFactors = FALSE)

Neetha.freq.1 <- Neetha.freq[which(Neetha.freq$`Number of fusions` < 250),]
Neetha.freq.1$quantile <- ifelse(Neetha.freq.1$`Number of fusions` <= quantile(Neetha.freq.1$`Number of fusions`, c(0.25)), "Low", ifelse(Neetha.freq.1$`Number of fusions` >=quantile(Neetha.freq.1$`Number of fusions`, c(0.75)), "Higher", "Middle"))
Neetha.freq.1.cBioPortal <- merge(Neetha.freq.1[,c("Number of fusions", "quantile")], BRCA.Clinical.cBioPortal, by="row.names")

### Fisher's exact test for HER2 positive
## By using only recurrent fusion data
n <- ncol(Neetha.freq.1.cBioPortal)
m <- nrow(Neetha.freq.1.cBioPortal[Neetha.freq.1.cBioPortal$IHC_HER2=="Positive",])
p <- nrow(Neetha.freq.1.cBioPortal[Neetha.freq.1.cBioPortal$quantile=="Higher",])
a <- nrow(Neetha.freq.1.cBioPortal[Neetha.freq.1.cBioPortal$IHC_HER2=="Positive" & Neetha.freq.1.cBioPortal$quantile=="Higher",])
b <- m-a
c <- p-a
d <- n-(a+b+c)
FisherTest <- matrix(c(a,b,c,d), nrow = 2, dimnames = list(c("Recurrent high fusion","Not high (recurrent fusion)"),c("HER2+ (recurrent fusion)", "Not-HER2+ (recurrent fusion)")))
fisher.test(FisherTest)
chisq.test(FisherTest)
Barnard::barnard.test(a,b,c,d)


### Fisher's exact test for HER2 positive
## By using recurrent fusion data and whole TCGA clinical data as suggested in IPA and other pathway tools
# n <- nrow(BRCA.Clinical.cBioPortal)
# m <- nrow(Neetha.freq.1.cBioPortal)
# p <- nrow(BRCA.Clinical.cBioPortal[BRCA.Clinical.cBioPortal$IHC_HER2=="Positive",])
# a <- nrow(Neetha.freq.1.cBioPortal[Neetha.freq.1.cBioPortal$IHC_HER2=="Positive" & Neetha.freq.1.cBioPortal$quantile=="Higher",])
# b <- p-a
# c <- m-a
# d <- n-(a+b+c)
# FisherTest <- matrix(c(a,b,c,d), nrow = 2, dimnames = list(c("Recurrent high fusion","Not recurrent high fusion"),c("HER2+ (TCGA)", "Not-HER2+ (TCGA)")))
# fisher.test(FisherTest)
# chisq.test(FisherTest)
# Barnard::barnard.test(a,b,c,d)