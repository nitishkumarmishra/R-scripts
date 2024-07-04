setwd("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/GDC-Harmonized/BRCA")
load("Comparative_Fusion.RData")
rownames(BRCA.Clinical) <- BRCA.Clinical$submitter_id
BRCA.Clinical.cBioPortal.select <- subset(BRCA.Clinical.cBioPortal, select = c(ER_STATUS_BY_IHC, PR_STATUS_BY_IHC, IHC_HER2, MENOPAUSE_STATUS))
BRCA.Clinical.cBioPortal.select$MENOPAUSE_STATUS <- gsub("[^0-9A-Za-z///' ]","" , BRCA.Clinical.cBioPortal.select$MENOPAUSE_STATUS ,ignore.case = TRUE)
BRCA.Clinical.cBioPortal.select$MENOPAUSE_STATUS <- sapply(strsplit(BRCA.Clinical.cBioPortal.select$MENOPAUSE_STATUS, " "),'[[', 1)
BRCA.Clinical.cBioPortal.select$MENOPAUSE_STATUS <- gsub("Not", "Not applicable", BRCA.Clinical.cBioPortal.select$MENOPAUSE_STATUS)
BRCA.Clinical.cBioPortal.select <- merge(BRCA.Clinical.cBioPortal.select, BRCA.Clinical, by = "row.names")
rownames(BRCA.Clinical.cBioPortal.select) <- BRCA.Clinical.cBioPortal.select$Row.names; BRCA.Clinical.cBioPortal.select$Row.names <- NULL
Neetha.freq.1 <- Neetha.freq[which(Neetha.freq$`Number of fusions` < 250),]
Neetha.freq.1$quantile <- ifelse(Neetha.freq.1$`Number of fusions` <= quantile(Neetha.freq.1$`Number of fusions`, c(0.25)), "Low", ifelse(Neetha.freq.1$`Number of fusions` >=quantile(Neetha.freq.1$`Number of fusions`, c(0.75)), "High", "Middle"))
Neetha.freq.1.cBioPortal <- merge(Neetha.freq.1[,c("Number of fusions", "quantile")], BRCA.Clinical.cBioPortal.select, by="row.names")
rownames(Neetha.freq.1.cBioPortal) <- Neetha.freq.1.cBioPortal$Row.names; Neetha.freq.1.cBioPortal$Row.names <- NULL
########### Fisher's exact test, by using only recurrent fusion data #########
## Clinical fewatures: ER_STATUS_BY_IHC, PR_STATUS_BY_IHC, IHC_HER2
## MolecularType= "Positive", "Negative"
## recurrentFusion= "High", "Low", "Middle"
Fisher_test <- function(ClinicalFeature, MolecularType, recurrentFusion){ 
  ### n:: total samples with recurrent fusion
  ### m:: recurrent fusion sample with HER2+/ER+/PR+ status
  ### p:: recurrent sample with "Higher" frequency
  ### a:: recurrent sample with "Higher" frequency and ER+/HER2+/PR+ status
  n <- nrow(Neetha.freq.1.cBioPortal) 
  m <- sum(Neetha.freq.1.cBioPortal[,ClinicalFeature]==MolecularType) 
  p <- nrow(Neetha.freq.1.cBioPortal[Neetha.freq.1.cBioPortal$quantile==recurrentFusion,]) 
  a <- sum(Neetha.freq.1.cBioPortal[,ClinicalFeature]==MolecularType&Neetha.freq.1.cBioPortal$quantile==recurrentFusion)
  b <- p-a; c <- m-a; d <- n-(a+b+c)
  FisherTest <- matrix(c(a,b,c,d), nrow = 2, byrow = TRUE, dimnames = list(c("recurrentFusion","Not-recurrentFusion"),c("ClinicalFeatureMolecularType", "Not-ClinicalFeatureMolecularType")))
  f1 <- fisher.test(FisherTest) #c1<- chisq.test(FisherTest)
  b1 <- Barnard::barnard.test(a,b,c,d)## In Barnard, barnard.test use ncol=2, while we are using nrow=2.
  ret <- data.frame(Fisher.P.value = f1$p.value, Barnard.one.p=b1$p.value[1], Barnard.two.p=b1$p.value[2])
  #return(list(contingency=FisherTest, Fisher.p.value=f1$p.value, Barnard.p.val.one.sided=b1$p.value[1], Barnard.p.val.two.sided=p.value[2]))
  return(list(contingency=FisherTest, Fisher.p.value=f1$p.value, Barnard.p.val.one=b1$p.value[1], Barnard.p.val.two=b1$p.value[2]))
  #return(ret)
  }
########## Command to run function ##########################
Fisher <- Fisher_test( recurrentFusion="High", ClinicalFeature="IHC_HER2", MolecularType="Positive")
Fisher <- Fisher_test( recurrentFusion="Low", ClinicalFeature="IHC_HER2", MolecularType="Positive")

Fisher <- Fisher_test( recurrentFusion="High", ClinicalFeature="ER_STATUS_BY_IHC", MolecularType="Positive")
Fisher <- Fisher_test( recurrentFusion="Low", ClinicalFeature="ER_STATUS_BY_IHC", MolecularType="Positive")

Fisher <- Fisher_test( recurrentFusion="High", ClinicalFeature="PR_STATUS_BY_IHC", MolecularType="Positive") ## It's p-value is OK
Fisher <- Fisher_test( recurrentFusion="Low", ClinicalFeature="PR_STATUS_BY_IHC", MolecularType="Positive")

Fisher <- Fisher_test( recurrentFusion="High", ClinicalFeature="MENOPAUSE_STATUS", MolecularType="Post")
Fisher <- Fisher_test( recurrentFusion="Low", ClinicalFeature="MENOPAUSE_STATUS", MolecularType="Post")

Fisher <- Fisher_test( recurrentFusion="High", ClinicalFeature="MENOPAUSE_STATUS", MolecularType="Pre")
Fisher <- Fisher_test( recurrentFusion="Low", ClinicalFeature="MENOPAUSE_STATUS", MolecularType="Pre")

######################## Save data #########################
save.image(file = "BRCA_Fishers.RData")
