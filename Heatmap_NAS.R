setwd("D:/OneDrive - University of Nebraska Medical Center/NAS")
load("ChAMP_DMP_DMR.RData")
rm(list=setdiff(ls(), c("myDMP", "myDMP.1", "myDMP.2",  "myNorm.PBC")))
################################################################
#https://jbengler.github.io/tidyheatmap/articles/tidyheatmap.html
library(tidyheatmap); library(dplyr); library(tidyverse); library(matrixStats)

############################################
pd <- ifelse(grepl("noOPIOD", colnames(myNorm.PBC)), "Group3", ifelse(grepl("OPIOD_No_NAS", colnames(myNorm.PBC)), "Group2", "Group1"))

Group2_to_Group1_beta_0.1 <- myDMP$Group1_to_Group2
#Group3_to_Group2_beta_0.1 <- myDMP$Group3_to_Group2[abs(myDMP$Group3_to_Group2$deltaBeta) >= 0.01,]
#Group3_to_Group2_beta_0.1 <- Group3_to_Group2_beta_0.1[order(abs(Group3_to_Group2_beta_0.1$deltaBeta), decreasing = TRUE),]
dataset1 <- myNorm.PBC[rownames(Group2_to_Group1_beta_0.1), ]
dataset1 <- dataset1[,grep("Group3", pd, invert = TRUE)]

dataset1 <- merge(dataset1, Group2_to_Group1_beta_0.1, by="row.names")
dataset1$Status <- ifelse(dataset1$deltaBeta >0 , "Hyper", "Hypo")
colnames(dataset1) <- gsub("OPIOD_", "", colnames(dataset1))
group <- pd[grep("Group3", pd, invert = TRUE)]
####################################
Num_Gene <- 120
Exp <- dataset1 %>%
  slice_max(abs(deltaBeta), n=Num_Gene) %>%
  mutate(log2BH= (log2BH = -log2(adj.P.Val))) %>%
  mutate(Feature = feature, "Delta Beta"=deltaBeta, "log2(BH)"=log2BH, CGI=cgi)%>%
  pivot_longer(grep("NAS", names(dataset1)),names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(group, Num_Gene))#rep(condition, 7); its seven gene so seventies of condition

ann_colors <- list(Group = c(Group1 = "red", Group2 = "blue"),
                   Status=c(Hyper="red", Hypo="blue"),
                   CGI=c(island="blue", opensea="red", shelf="green", shore="darkcyan"),
                   logFC= c("blue","gray","red"),
                   "log2(BH)"=c("blue","gray","red"),
                   "Delta Beta"=c("blue","gray","red"),
                   Feature=c(Body="darkgreen",IGR="coral2",ExonBnd="cyan4",TSS1500="red","5'UTR"="orange","3'UTR"= "magenta4",TSS200="blue","1stExon"="blueviolet"))

tidy_heatmap(Exp,
             rows = "Row.names",
             columns = Sample,
             values = Expression,
             scale = "none",
             annotation_col = c(Group),
             annotation_row = c(Status, "Delta Beta",logFC, "log2(BH)",CGI, Feature),
             cluster_cols = TRUE,
             clustering_method = "ward.D2",
             #clustering_distance_cols = "euclidean",
             clustering_distance_cols = "manhattan",
             cluster_rows = TRUE,
             color_legend_n = 7,
             colors = c("red","gray","green"),
             annotation_colors = ann_colors,
             fontsize = 10,
             fontsize_row = 6,
             fontsize_col = 8,
             #angle_col = 45,
             height = 12,
             width = 10,
             #filename = "Heatmap_grou1_vs2_120probes.pdf"
)

####################################################################
####################################################################
pd.2 <- gsub("Group2", "Group1", pd)
Group1_and_2_vs_Group3 <- myDMP.2$Group3_to_Group1

dataset1 <- myNorm.PBC[rownames(Group1_and_2_vs_Group3), ]
#dataset1 <- dataset1[,grep("Group3", pd, invert = TRUE)]

dataset1 <- merge(dataset1, Group1_and_2_vs_Group3, by="row.names")
dataset1$Status <- ifelse(dataset1$deltaBeta >0 , "Hyper", "Hypo")
colnames(dataset1) <- gsub("OPIOD_", "", colnames(dataset1))
#group <- pd[grep("Group3", pd, invert = TRUE)]
group <- pd

Num_Gene <- 116
Exp <- dataset1 %>%
  slice_max(abs(deltaBeta), n=Num_Gene) %>%
  mutate(log2BH= (log2BH = -log2(adj.P.Val))) %>%
  mutate(Feature = feature, "Delta Beta"=deltaBeta, "log2(BH)"=log2BH, CGI=cgi)%>%
  pivot_longer(grep("NAS", names(dataset1)),names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(group, Num_Gene))#rep(condition, 7); its seven gene so seventies of condition

ann_colors <- list(Group = c(Group1 = "red", Group2 = "blue", Group3="green"),
                   Status=c(Hyper="red", Hypo="blue"),
                   CGI=c(island="blue", opensea="red", shelf="green", shore="darkcyan"),
                   logFC= c("blue","gray","red"),
                   "log2(BH)"=c("blue","gray","red"),
                   "Delta Beta"=c("blue","gray","red"),
                   Feature=c(Body="darkgreen",IGR="coral2",ExonBnd="cyan4",TSS1500="red","5'UTR"="orange","3'UTR"= "magenta4",TSS200="blue","1stExon"="blueviolet"))

tidy_heatmap(Exp,
             rows = "Row.names",
             columns = Sample,
             values = Expression,
             scale = "none",
             annotation_col = c(Group),
             annotation_row = c(Status, "Delta Beta",logFC, "log2(BH)",CGI, Feature),
             cluster_cols = TRUE,
             clustering_method = "ward.D2",
             clustering_distance_cols = "euclidean",
             #clustering_distance_cols = "manhattan",
             cluster_rows = TRUE,
             color_legend_n = 7,
             colors = c("red4","red1","gray","green1", "green4"),
             annotation_colors = ann_colors,
             fontsize = 10,
             fontsize_row = 6,
             fontsize_col = 5,
             #angle_col = 45,
             height = 13,
             width = 11,
             filename = "Heatmap1_and_2_vs_Group3_probe116.pdf"
)


####################################################################
####################################################################
# Group1 vs Group3
pd.3 <- grep("Group2", pd, invert = TRUE)
Group3_vs_Group1 <- myDMP$Group3_to_Group1

dataset1 <- myNorm.PBC[rownames(Group3_vs_Group1), pd.3]
#dataset1 <- dataset1[,grep("Group3", pd, invert = TRUE)]

dataset1 <- merge(dataset1, Group3_vs_Group1, by="row.names")
dataset1$Status <- ifelse(dataset1$deltaBeta >0 , "Hyper", "Hypo")
colnames(dataset1) <- gsub("OPIOD_", "", colnames(dataset1))
#group <- pd[grep("Group3", pd, invert = TRUE)]
group <- pd[pd.3]

Num_Gene <- 290
Exp <- dataset1 %>%
  slice_max(abs(deltaBeta), n=Num_Gene) %>%
  mutate(log2BH= (log2BH = -log2(adj.P.Val))) %>%
  mutate(Feature = feature, "Delta Beta"=deltaBeta, "log2(BH)"=log2BH, CGI=cgi)%>%
  pivot_longer(grep("NAS", names(dataset1)),names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(group, Num_Gene))#rep(condition, 7); its seven gene so seventies of condition

ann_colors <- list(Group = c(Group1 = "red", Group3="green"),
                   Status=c(Hyper="red", Hypo="blue"),
                   CGI=c(island="blue", opensea="red", shelf="green", shore="darkcyan"),
                   logFC= c("blue","gray","red"),
                   "log2(BH)"=c("blue","gray","red"),
                   "Delta Beta"=c("blue","gray","red"),
                   Feature=c(Body="darkgreen",IGR="coral2",ExonBnd="cyan4",TSS1500="red","5'UTR"="orange","3'UTR"= "magenta4",TSS200="blue","1stExon"="blueviolet"))

tidy_heatmap(Exp,
             rows = "Row.names",
             columns = Sample,
             values = Expression,
             scale = "none",
             annotation_col = c(Group),
             annotation_row = c(Status, "Delta Beta",logFC, "log2(BH)",CGI, Feature),
             cluster_cols = TRUE,
             clustering_method = "ward.D2",
             clustering_distance_cols = "euclidean",
             #clustering_distance_cols = "manhattan",
             cluster_rows = TRUE,
             color_legend_n = 7,
             colors = c("red4","red1","gray","green1", "green4"),
             annotation_colors = ann_colors,
             fontsize = 10,
             fontsize_row = 6,
             show_rownames = FALSE,
             fontsize_col = 5,
             #angle_col = 45,
             height = 13,
             width = 11,
             filename = "Heatmap1_vs_Group3_probe290.pdf"
)
#####################################################################
save.image("NAS_Heatmap.RData")
