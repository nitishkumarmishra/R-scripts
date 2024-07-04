library(TCGAbiolinks)
library(dplyr)
library(DT)
# ECD (ENSG00000122882) ,ERBB2 (ENSG00000141736), ESR1 (ENSG00000091831), ESR2 (ENSG00000140009), PGR (ENSG00000082175)
setwd("D:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/GDC-Harmonized/BRCA/")
load("BRCA.RData")
table(substr(colnames(BRCA.FPKM), 14, 16))
#  01A  01B  01C  06A  11A  11B 
# 1077   24    1    7   99   14
FPKM_Tumor <- BRCA.FPKM %>% 
  dplyr::select(grep("-01", names(BRCA.FPKM)))
## Remove duplicated samples
colnames(FPKM_Tumor) <- substr(colnames(FPKM_Tumor), 1, 15)
FPKM_Tumor <- FPKM_Tumor[,!duplicated(colnames(FPKM_Tumor))]

FPKM_Normal <- BRCA.FPKM %>% 
  dplyr::select(grep("-11A|-11B", names(BRCA.FPKM)))
## Remove duplicated samples
colnames(FPKM_Normal) <- substr(colnames(FPKM_Normal), 1, 15)
FPKM_Normal <- FPKM_Normal[,!duplicated(colnames(FPKM_Normal))]

FPKM_Metastasis <- BRCA.FPKM %>% 
  dplyr::select(grep("-06A", names(BRCA.FPKM)))
colnames(FPKM_Metastasis) <- substr(colnames(FPKM_Metastasis), 1, 15)
FPKM_Metastasis <- FPKM_Metastasis[,!duplicated(colnames(FPKM_Metastasis))]
#impute::impute.knn() does't work, because of memory size
#FPKM_Normal.knn <- impute::impute.knn(as.matrix(FPKM_Normal), k = 15, rowmax = 0.8, colmax = 0.8, maxp = 500, rng.seed=123)
rowMeans(FPKM_Normal[grep("ENSG00000122882", rownames(FPKM_Normal)),], na.rm = TRUE)
# 10.428
ECD <- FPKM_Tumor[grep("ENSG00000122882", rownames(FPKM_Tumor)),]
ECD.Status <- ifelse(ECD >10.428, "ECD_High", "ECD_Low")
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

Surviva_Cell.2018 <- readr::read_delim("TCGA Survival Cell 2018.txt", delim = "\t")%>%
  filter(type == "BRCA") %>%
  select(bcr_patient_barcode, starts_with("OS"))

Survival_Data <- inner_join(cBioPortal.clinical, Surviva_Cell.2018, by ="bcr_patient_barcode")
Survival_Data <- Survival_Data[-grep("#N/A", Survival_Data$OS.time),]

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
  rename(ECD = ENSG00000122882.9,ERBB2 = ENSG00000141736.12, ESR1 = ENSG00000091831.20, ESR2 = ENSG00000140009.17, PGR = ENSG00000082175.13)
#################### correlation analysis ###############
library(ggpubr)
ECD.Exp.log2 <- ECD.Exp %>%
  mutate(ECD = log2(ECD.Exp$ECD+1)) %>%
  mutate(ERBB2 = log2(ECD.Exp$ERBB2+1)) %>%
  mutate(ESR1 = log2(ECD.Exp$ESR1+1)) %>%
  mutate(ESR2 = log2(ECD.Exp$ESR2+1)) %>%
  mutate(PGR = log2(ECD.Exp$PGR+1))
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
## Spearnan's correlation analysis ##
ggscatter(ECD.Exp.log2, x = "ECD", y = "ERBB2",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman", label.x = 3, label.y = 30)  # Add correlation coefficient
dev.print(pdf, 'ECD and ERBB2 Spearman correlation.pdf', width = 8, height = 8)

###### t test to check difference ##########
###### It doesn't make sense, both gene have different level of expression
ECD.Exp.log2.melt <- reshape2::melt(ECD.Exp.log2) %>%
  rename(Gene=variable, Expression=value) %>%
  filter(grepl('ECD|ERBB2', Gene))

p <- ggboxplot(ECD.Exp.log2.melt, x = "Gene", y = "Expression", color = "Gene", palette = "jco", add = "jitter")
p + stat_compare_means(method = "t.test")
dev.print(pdf, 'ECD and ERBB2 t-test.pdf', width = 8, height = 8)
##########################################################
##########################################################
# Survival analysis ############
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)

Survival_Data %>%
  analyse_survival(vars(OS.time, OS), by=ECD_Status) ->  result

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
dev.print(pdf, 'ECD Survival plot.pdf', width = 10, height = 8)

##########################################################
####### Patients with over 100 days survival time ########
#Survival_Data %>%
  #filter(OS.time > 100) %>%count (PR_Status) # Table of each PR group


Survival_Data.select <- Survival_Data %>%
  filter(OS.time > 100) %>%
  analyse_survival(vars(OS.time, OS), by=ECD_Status) ->  result

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
dev.print(pdf, 'ECD OS Over 100 Days Survival plot.pdf', width = 10, height = 8)
##########################################################
############# Survival analysis in TNBC ##################
Survival_Data %>%
  filter(OS.time > 100) %>%
  filter(TNBC=="TNBC") %>%
  analyse_survival(vars(OS.time, OS), by=ECD_Status) ->  result
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
dev.print(pdf, 'ECD OS Over 100 Days in TNBC Survival plot.pdf', width = 10, height = 8)
############# Survival analysis in ER+/Her2+ ################
Survival_Data %>%
  filter(OS.time > 100) %>%
  filter(ER_Status=="Positive") %>%
  analyse_survival(vars(OS.time, OS), by=ECD_Status) ->  result
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
dev.print(pdf, 'ECD OS Over 100 Days in ER+ Survival plot.pdf', width = 10, height = 8)

############# Survival analysis in PR+ #################
Survival_Data %>%
  filter(OS.time > 100) %>%
  filter(PR_Status=="Positive") %>%
  analyse_survival(vars(OS.time, OS), by=ECD_Status) ->  result
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
dev.print(pdf, 'ECD OS Over 100 Days in PR+ Survival plot.pdf', width = 10, height = 8)
##########################################################
######## Correlation between ECD and ERBB2 in TNBC #######
TNBC_ECD.exp <- 
  inner_join(Survival_Data, ECD.Exp.log2, by = "bcr_patient_barcode") %>%
  filter(TNBC=="TNBC")

ggscatter(TNBC_ECD.exp, x = "ECD", y = "ERBB2",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman", label.x = 3, label.y = 30)  # Add correlation coefficient
dev.print(pdf, 'ECD and ERBB2 in TNBC Spearman correlation.pdf', width = 8, height = 8)

##########################################################
######## Correlation between ECD and PR in TNBC #######
TNBC_ECD.exp <- 
  inner_join(Survival_Data, ECD.Exp.log2, by = "bcr_patient_barcode") %>%
  filter(TNBC=="TNBC")

ggscatter(TNBC_ECD.exp, x = "ECD", y = "ERBB2",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman", label.x = 3, label.y = 30)  # Add correlation coefficient
dev.print(pdf, 'ECD and ERBB2 in TNBC Spearman correlation.pdf', width = 8, height = 8)


##########################################################
save.image("ECD_Band.RData")
##########################################################
# Sample with more than one tumor data (01A/01B/01C)
head(BRCA.FPKM %>%
  dplyr::select(grep("TCGA-A7-A13E", names(BRCA.FPKM))))
head(BRCA.FPKM %>%
       dplyr::select(grep("TCGA-A7-A0DB", names(BRCA.FPKM))))
###########################################################
###########################################################
mean(FPKM_Normal[grep("ENSG00000122882", rownames(FPKM_Normal)),])
