library(TCGAbiolinks)
library(dplyr)
library(DT)
# ECD (ENSG00000122882) ,ERBB2 (ENSG00000141736), ESR1 (ENSG00000091831), ESR2 (ENSG00000140009), PGR (ENSG00000082175)
setwd("D:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/GDC-Harmonized/BRCA/")
load("BRCA.RData")
table(substr(colnames(BRCA.FPKM.UQ), 14, 16))
#  01A  01B  01C  06A  11A  11B 
# 1077   24    1    7   99   14
FPKM_Tumor <- BRCA.FPKM.UQ %>% 
  dplyr::select(grep("-01", names(BRCA.FPKM.UQ)))
## Remove duplicated samples
colnames(FPKM_Tumor) <- substr(colnames(FPKM_Tumor), 1, 15)
FPKM_Tumor <- FPKM_Tumor[,!duplicated(colnames(FPKM_Tumor))]

FPKM_Normal <- BRCA.FPKM.UQ %>% 
  dplyr::select(grep("-11A|-11B", names(BRCA.FPKM.UQ)))
## Remove duplicated samples
colnames(FPKM_Normal) <- substr(colnames(FPKM_Normal), 1, 15)
FPKM_Normal <- FPKM_Normal[,!duplicated(colnames(FPKM_Normal))]

FPKM_Metastasis <- BRCA.FPKM.UQ %>% 
  dplyr::select(grep("-06A", names(BRCA.FPKM.UQ)))
colnames(FPKM_Metastasis) <- substr(colnames(FPKM_Metastasis), 1, 15)
FPKM_Metastasis <- FPKM_Metastasis[,!duplicated(colnames(FPKM_Metastasis))]
#impute::impute.knn() does't work, because of memory size
#FPKM_Normal.knn <- impute::impute.knn(as.matrix(FPKM_Normal), k = 15, rowmax = 0.8, colmax = 0.8, maxp = 500, rng.seed=123)
rowMeans(FPKM_Normal[grep("ENSG00000122882", rownames(FPKM_Normal)),], na.rm = TRUE)
# 10.428
ECD <- FPKM_Tumor[grep("ENSG00000122882", rownames(FPKM_Tumor)),]
ECD.Status <- ifelse(ECD >median(as.matrix(ECD)), "ECD_High", "ECD_Low")
ECD.High <- grep("ECD_High", ECD.Status)
ECD.High.TCGA_ID <- names(ECD[ECD.High])
ECD.Low <- grep("ECD_Low", ECD.Status)
ECD.Low.TCGA_ID <- names(ECD[ECD.Low])


##########################################################
##########################################################
cBioPortal.clinical <- readr::read_delim("clinical_PANCAN_patient_with_followup.tsv", delim = "\t") %>%
  filter(acronym == "BRCA") %>%
  select(bcr_patient_barcode, gender, vital_status, radiation_therapy, days_to_death, days_to_last_followup,age_at_initial_pathologic_diagnosis, race, ethnicity, ends_with("receptor_status"), -starts_with("metastatic_breast")) %>%
  rename(PR_Status = breast_carcinoma_progesterone_receptor_status, ER_Status= breast_carcinoma_estrogen_receptor_status, HER2_Status = lab_proc_her2_neu_immunohistochemistry_receptor_status) %>%
  mutate(TNBC = if_else(PR_Status=="Negative" & ER_Status == "Negative" & HER2_Status=="Negative", "TNBC", "NO")) %>%
  mutate(Time=if_else(vital_status =="Alive", days_to_last_followup, days_to_death))
#mutate(Time=as.numeric(Time)*0.032854884083862)

#tmp <- readr::read_delim("brca_tcga_pan_can_atlas_2018.tar/data_clinical_patient.txt", delim = "\t") %>%
#  select(PATIENT_ID,	SUBTYPE,CANCER_TYPE_ACRONYM, starts_with("OS_"))

Survival_Cell.2018 <- readr::read_delim("TCGA Survival Cell 2018.txt", delim = "\t")%>%
  filter(type == "BRCA") %>%
  select(bcr_patient_barcode, starts_with("OS"))

Survival_Data <- inner_join(cBioPortal.clinical, Survival_Cell.2018, by ="bcr_patient_barcode")
Survival_Data <- Survival_Data[-grep("#N/A", Survival_Data$OS.time),]


Survival_PANCA.2018 <- readr::read_delim("brca_tcga_pan_can_atlas_2018_clinical_data.tsv", delim = "\t")%>%
  select("Patient ID", Subtype) %>%
  rename(bcr_patient_barcode="Patient ID")

Survival_Data  <- left_join(Survival_Data, Survival_PANCA.2018, by ="bcr_patient_barcode")

Survival_Data  <- Survival_Data %>%
  mutate(ECD_Status = if_else(bcr_patient_barcode %in% substr(ECD.High.TCGA_ID, 1, 12), "High", "Low")) %>%
  mutate(OS.time=as.numeric(OS.time))%>%
  mutate(time=as.numeric(OS.time)*0.032854884083862)
## Remove TCGA ID from Clnical data which don't have RNAseq data
Survival_Data <- Survival_Data[Survival_Data$bcr_patient_barcode %in% substr(c(ECD.High.TCGA_ID, ECD.Low.TCGA_ID),1, 12),]
colnames(FPKM_Tumor) <- substr(colnames(FPKM_Tumor),1, 12)
FPKM_Tumor <- FPKM_Tumor[,Survival_Data$bcr_patient_barcode]
FPKM_Tumor <- t(FPKM_Tumor)
#head(FPKM_Tumor[,grep("ENSG00000122882|ENSG00000141736|ENSG00000091831|ENSG00000140009|ENSG00000082175", colnames(FPKM_Tumor))])
ECD.Exp <- as_tibble(FPKM_Tumor, rownames="rownames") %>%
  select(bcr_patient_barcode=rownames, ENSG00000122882.9, ENSG00000141736.12, ENSG00000091831.20, ENSG00000140009.17, ENSG00000082175.13)

ECD.Exp <- ECD.Exp %>%
  rename(ECD = ENSG00000122882.9, ERBB2 = ENSG00000141736.12, ESR1 = ENSG00000091831.20, ESR2 = ENSG00000140009.17, PGR = ENSG00000082175.13)


# ECD.Exp.N <- as_tibble(FPKM_Normal, rownames="rownames") %>%
#   select(bcr_patient_barcode=rownames, ENSG00000122882.9, ENSG00000141736.12, ENSG00000091831.20, ENSG00000140009.17, ENSG00000082175.13)
# ECD.Exp.N <- ECD.Exp.N %>%
#   rename(ECD = ENSG00000122882.9,ERBB2 = ENSG00000141736.12, ESR1 = ENSG00000091831.20, ESR2 = ENSG00000140009.17, PGR = ENSG00000082175.13)


####################################################
colnames(FPKM_Normal) <- substr(colnames(FPKM_Normal),1, 12)
FPKM_Normal <- FPKM_Normal[ , colnames(FPKM_Normal) %in% Survival_Data$bcr_patient_barcode]
FPKM_Normal <- t(FPKM_Normal)

ECD.Exp.N <- as_tibble(FPKM_Normal, rownames="rownames") %>%
  select(bcr_patient_barcode=rownames, ENSG00000122882.9, ENSG00000141736.12, ENSG00000091831.20, ENSG00000140009.17, ENSG00000082175.13)
ECD.Exp.N <- ECD.Exp.N %>%
  rename(ECD = ENSG00000122882.9,ERBB2 = ENSG00000141736.12, ESR1 = ENSG00000091831.20, ESR2 = ENSG00000140009.17, PGR = ENSG00000082175.13)

#################### log2 transformation ###############
ECD.Exp.log2 <- ECD.Exp %>%
  mutate(ECD = log2(ECD.Exp$ECD+1)) %>%
  mutate(ERBB2 = log2(ECD.Exp$ERBB2+1)) %>%
  mutate(ESR1 = log2(ECD.Exp$ESR1+1)) %>%
  mutate(ESR2 = log2(ECD.Exp$ESR2+1)) %>%
  mutate(PGR = log2(ECD.Exp$PGR+1))

ECD.Exp.N.log2 <- ECD.Exp.N %>%
  mutate(ECD = log2(ECD.Exp.N$ECD+1)) %>%
  mutate(ERBB2 = log2(ECD.Exp.N$ERBB2+1)) %>%
  mutate(ESR1 = log2(ECD.Exp.N$ESR1+1)) %>%
  mutate(ESR2 = log2(ECD.Exp.N$ESR2+1)) %>%
  mutate(PGR = log2(ECD.Exp.N$PGR+1))

ECD.Exp.N.log2 <- ECD.Exp.N.log2 %>%
  mutate(Status=rep("Normal"))
ECD.Exp.log2 <- ECD.Exp.log2 %>%
  mutate(Status=rep("Tumor"))

ECD.data <- rbind(ECD.Exp.N.log2, ECD.Exp.log2)
#################### correlation analysis ###############
library(ggpubr)
#################### log2 transformation ###############
#########################################################
## Pearson's correlation analysis ##
ggscatter(ECD.Exp.log2, x = "ECD", y = "ERBB2",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "pearson", label.x = 3, label.y = 30)  # Add correlation coefficient
dev.print(pdf, 'ECD and ERBB2 Pearson correlation.pdf', width = 8, height = 8)

#########################################################
####### ECD HER2 correlation in HER2+ (ERBB2+) ##########
#########################################################
Survival_Data.HER2 <- Survival_Data %>%
  filter(HER2_Status =="Positive")
ECD.Exp.log2.HER2 <- ECD.Exp.log2[ECD.Exp.log2$bcr_patient_barcode %in% Survival_Data.HER2$bcr_patient_barcode, ]
#ECD.data.HER2 <- rbind(ECD.Exp.N.log2, ECD.Exp.log2.HER2)

ggscatter(ECD.Exp.log2.HER2, x = "ECD", y = "ERBB2",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "pearson", label.x = 3, label.y = 30)  # Add correlation coefficient
dev.print(pdf, 'ECD and ERBB2 Pearson correlation in HER2+.pdf', width = 8, height = 8)
#########################################################
#########################################################
################ ECD tumor normal t-test ################
#########################################################

ggboxplot(ECD.data, x = "Status", y = "ECD",
          title = "ECD", ylab = "Expression",
          color = "Status", palette = "lancet",
          add = "jitter",                              # Add jittered points
          add.params = list(size = 0.1, jitter = 0.2))+ 
  stat_compare_means(method = "t.test")
#stat_compare_means(method = "wilcox.test")
#stat_compare_means(method = "anova")

dev.print(pdf, 'ECD Tumor vs Normal t-test.pdf', width = 8, height = 8)
#########################################################
################## ECD ER+/PR+ t-test ###################
#########################################################
Survival_Data <- Survival_Data %>%
  mutate(ErPr_Status=ifelse(PR_Status=="Positive" & ER_Status == "Positive", "Positive", "Negative"))
Survival_Data.ErPr <- Survival_Data %>%
  filter(ErPr_Status =="Positive")
ECD.Exp.log2.ErPr <- ECD.Exp.log2[ECD.Exp.log2$bcr_patient_barcode %in% Survival_Data.ErPr$bcr_patient_barcode, ]
ECD.data.ErPr <- rbind(ECD.Exp.N.log2, ECD.Exp.log2.ErPr)

ggboxplot(ECD.data.ErPr, x = "Status", y = "ECD",
          title = "ER+/PR+", ylab = "Expression",
          color = "Status", palette = "lancet",
          add = "jitter",                              # Add jittered points
          add.params = list(size = 0.1, jitter = 0.2))+ 
  stat_compare_means(method = "t.test")
#stat_compare_means(method = "wilcox.test")
#stat_compare_means(method = "anova")

dev.print(pdf, 'ECD ER+PR+ vs Normal t-test.pdf', width = 8, height = 8)

#########################################################
############## ECD HER2+ (ERBB2+) t-test ################
#########################################################
Survival_Data.HER2 <- Survival_Data %>%
  filter(HER2_Status =="Positive")
ECD.Exp.log2.HER2 <- ECD.Exp.log2[ECD.Exp.log2$bcr_patient_barcode %in% Survival_Data.HER2$bcr_patient_barcode, ]
ECD.data.HER2 <- rbind(ECD.Exp.N.log2, ECD.Exp.log2.HER2)

ggboxplot(ECD.data.HER2, x = "Status", y = "ECD",
          title = "HER2+ (ERBB2+)", ylab = "Expression",
          color = "Status", palette = "lancet",
          add = "jitter",                              # Add jittered points
          add.params = list(size = 0.1, jitter = 0.2))+ 
  stat_compare_means(method = "t.test")
#stat_compare_means(method = "wilcox.test")
#stat_compare_means(method = "anova")

dev.print(pdf, 'ECD HER2+ vs Normal t-test.pdf', width = 8, height = 8)

#########################################################
#################### ECD TNBC t-test ####################
#########################################################
Survival_Data.TNBC <- Survival_Data %>%
  filter(TNBC =="TNBC")
ECD.Exp.log2.TNBC <- ECD.Exp.log2[ECD.Exp.log2$bcr_patient_barcode %in% Survival_Data.TNBC$bcr_patient_barcode, ]
ECD.data.TNBC <- rbind(ECD.Exp.N.log2, ECD.Exp.log2.TNBC)

ggboxplot(ECD.data.TNBC, x = "Status", y = "ECD",
          title = "TNBC", ylab = "Expression",
          color = "Status", palette = "lancet",
          add = "jitter",                              # Add jittered points
          add.params = list(size = 0.1, jitter = 0.2))+ 
  stat_compare_means(method = "t.test")
#stat_compare_means(method = "wilcox.test")
#stat_compare_means(method = "anova")

dev.print(pdf, 'ECD TNBC vs Normal t-test.pdf', width = 8, height = 8)



#########################################################
#################### TNBC and basal t-test ####################
#########################################################
Survival_Data <- Survival_Data %>%
  mutate(TNBC_Basal=ifelse(TNBC=="TNBC" & Subtype == "BRCA_Basal", "Positive", "Negative"))
Data_Survival <- inner_join(Survival_Data, ECD.Exp.log2)

Survival_Data.TNBC_Basal <- Survival_Data %>%
  filter(TNBC_Basal =="Positive")

Survival_Data.TNBC_No_Basal <- Survival_Data %>%
  filter(TNBC=="TNBC" & TNBC_Basal =="Negative")
ECD.Exp.N.log2_TNBC_NoBasal <- ECD.Exp.log2[ECD.Exp.log2$bcr_patient_barcode %in% Survival_Data.TNBC_No_Basal$bcr_patient_barcode, ]
ECD.Exp.N.log2_TNBC_NoBasal$Status <- gsub("Tumor", "Non-Basal", ECD.Exp.N.log2_TNBC_NoBasal$Status)


##ECD.Exp.log2.TNBC_Basal <- ECD.Exp.log2[ECD.Exp.log2$bcr_patient_barcode %in% Survival_Data.TNBC_Basal$bcr_patient_barcode, ]
#ECD.Exp.log2.TNBC_Basal$Status <- gsub("Tumor", "Basal", ECD.Exp.log2.TNBC_Basal$Status)

ECD.Exp.N.log2.1 <- ECD.Exp.N.log2
ECD.Exp.N.log2.1$Status <- gsub("Tumor", "Normal", ECD.Exp.N.log2.1$Status)


ECD.data.TNBC_Basal <- rbind(ECD.Exp.N.log2.1, ECD.Exp.log2.TNBC_Basal)

#ECD.data.TNBC_Basal <- rbind(ECD.Exp.N.log2_TNBC_NoBasal, ECD.Exp.log2.TNBC_Basal)

ggboxplot(ECD.data.TNBC_Basal, x = "Status", y = "ECD",
          title = "TNBC basal and normal", ylab = "Expression",
          color = "Status", palette = "lancet",
          add = "jitter",                              # Add jittered points
          add.params = list(size = 0.1, jitter = 0.2))+ 
  stat_compare_means(method = "t.test")

dev.print(pdf, 'ECD TNBC Basal vs Non-basal t-test.pdf', width = 8, height = 8)
dev.print(pdf, 'ECD TNBC Basal vs normal t-test.pdf', width = 8, height = 8)
##########################################################
############ ECD Overall Survival analysis ###############
##########################################################
###### ECD Overall Survival analysis in all Tumor ########
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)

library(survminer)

res.cut <- surv_cutpoint(Data_Survival, time = "OS.time", event = "OS",
                         variables = c("ECD", "ERBB2"), minprop = 0.25)
summary(res.cut)
plot(res.cut, "ECD", palette = "npg")

res.cat <- surv_categorize(res.cut)
res.cat$OS.time <- res.cat$OS.time*0.032854884083862

library("survival")
fit <- survfit(Surv(OS.time, OS) ~ECD, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = FALSE, pval = TRUE, palette =  "npg", xlab="Survival Months")
dev.print(pdf, 'ECD Survival plot FPKM-UQ.pdf', width = 12, height = 10)
######################################
##########################################################
res.cut <- surv_cutpoint(Data_Survival, time = "OS.time", event = "OS",
                         variables = c("ECD", "ERBB2"), minprop = 0.25)
res.cat <- surv_categorize(res.cut)

res.cat %>%
  analyse_survival(vars(OS.time, OS), by=ECD) ->  result

kaplan_meier_plot(result,
                  break.time.by="breakByYear",
                  xlab=".OS.months",
                  legend.title="ECD Status",
                  hazard.ratio=T,
                  risk.table=TRUE,
                  table.layout="clean",
                  ggtheme = ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                                           legend.key = element_rect(fill = "white")))
dev.print(pdf, 'ECD Survival plot FPKM-UQ Tidyvers.pdf', width = 12, height = 10)


######## ECD Overall Survival analysis in ERBB2+ #########
Data_Survival_HER2 <- Data_Survival %>%
  filter(HER2_Status=="Positive")
res.cut <- surv_cutpoint(Data_Survival_HER2, time = "OS.time", event = "OS",
                         variables = c("ECD"), minprop = 0.10)
summary(res.cut)
plot(res.cut, "ECD", palette = "npg")

res.cat <- surv_categorize(res.cut)
res.cat$OS.time <- res.cat$OS.time*0.032854884083862

fit <- survfit(Surv(OS.time, OS) ~ECD, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = FALSE, pval = TRUE, palette =  "npg", xlab="Survival Months")
dev.print(pdf, 'ECD Survival plot in ERBB2+ HER2+ FPKM-UQ.pdf', width = 12, height = 10)


######## ECD Overall Survival analysis in ER+/PR+ #########
Data_Survival_ER_and_PR <- Data_Survival %>%
  filter(ErPr_Status=="Positive")
res.cut <- surv_cutpoint(Data_Survival_ER_and_PR, time = "OS.time", event = "OS",
                         variables = c("ECD"), minprop = 0.10)
summary(res.cut)
plot(res.cut, "ECD", palette = "npg")

res.cat <- surv_categorize(res.cut)
res.cat$OS.time <- res.cat$OS.time*0.032854884083862

fit <- survfit(Surv(OS.time, OS) ~ECD, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = FALSE, pval = TRUE, palette =  "npg", xlab="Survival Months")
dev.print(pdf, 'ECD Survival plot in ER+ PR+ FPKM-UQ.pdf', width = 12, height = 10)


######## ECD Overall Survival analysis in TNBC #########
Data_Survival_TNBC <- Data_Survival %>%
  filter(TNBC=="TNBC")
res.cut <- surv_cutpoint(Data_Survival_TNBC, time = "OS.time", event = "OS",
                         variables = c("ECD"), minprop = 0.10)
summary(res.cut)
plot(res.cut, "ECD", palette = "npg")

res.cat <- surv_categorize(res.cut)
res.cat$OS.time <- res.cat$OS.time*0.032854884083862

fit <- survfit(Surv(OS.time, OS) ~ECD, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = FALSE, pval = TRUE, palette =  "npg", xlab="Survival Months")
dev.print(pdf, 'ECD Survival plot in TNBC FPKM-UQ.pdf', width = 12, height = 10)


######## ECD Overall Survival analysis in TNBC #########

Data_Survival_TNBC_Basal <- Data_Survival %>%
  filter(TNBC=="TNBC" & Subtype=="BRCA_Basal")

res.cut <- surv_cutpoint(Data_Survival_TNBC_Basal, time = "OS.time", event = "OS",
                         variables = c("ECD"), minprop = 0.10)
summary(res.cut)
plot(res.cut, "ECD", palette = "npg")

res.cat <- surv_categorize(res.cut)
res.cat$OS.time <- res.cat$OS.time*0.032854884083862

fit <- survfit(Surv(OS.time, OS) ~ECD, data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = FALSE, pval = TRUE, palette =  "npg", xlab="Survival Months")
dev.print(pdf, 'ECD Survival plot in TNBC and Basal FPKM-UQ.pdf', width = 12, height = 10)


#########################################################################
#########################################################################
################## Survival based on median and quartile ################
Data_Survival_ER_and_PR <- Data_Survival %>%
  filter(ErPr_Status=="Positive") %>%
  mutate(plot=ifelse(ECD > median(ECD), "High", "Low"))

fit <- survfit(Surv(OS.time, OS) ~plot, data = Data_Survival_ER_and_PR)
ggsurvplot(fit, risk.table = TRUE, conf.int = FALSE, pval = TRUE, palette =  "npg", xlab="Survival Months")
dev.print(pdf, 'ECD Survival plot in ER+ PR+ FPKM-UQ Median.pdf', width = 12, height = 10)



Data_Survival_ER_and_PR <- Data_Survival %>%
  filter(ErPr_Status=="Positive") %>%
  mutate(plot=ntile(ECD, 3)) %>%
  filter(plot != "2")

fit <- survfit(Surv(OS.time, OS) ~plot, data = Data_Survival_ER_and_PR)
ggsurvplot(fit, risk.table = TRUE, conf.int = FALSE, pval = TRUE, palette =  "npg", xlab="Survival Months")
dev.print(pdf, 'ECD Survival plot in ER+ PR+ FPKM-UQ Upper and Lower quartile.pdf', width = 12, height = 10)

##########################################################  
save.image("ECD_Bhopal_FPKM_UQ.RData")
##########################################################
