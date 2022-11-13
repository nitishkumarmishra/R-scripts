
library(readxl)
library(data.table)
library(ggplot2)
library(genefilter)
library(dplyr)

setwd("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot")

# Load data
load("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/NAS_withXY.RData")


rm(list=ls()[! ls() %in% c('myNorm.PBC')])


#Group1 <- unique(gsub("\\..*","",colnames(sample_set)[grepl("OPIOD_NAS-Aff", colnames(sample_set))]))
#Group2 <- unique(gsub("\\..*","",colnames(sample_set)[grepl("OPIOD_No_NAS-Contr", colnames(sample_set))]))



# OPOID NAS
Group1 <- grepl("OPIOD_NAS-Aff", colnames(myNorm.PBC))
Group1 <- myNorm.PBC[,Group1]

# OPIOD_No_NAS
Group2 <- grepl("\\bOPIOD_No_NAS-Contr\\b", colnames(myNorm.PBC))
Group2 <- myNorm.PBC[,Group2]

#noOPOID
Group3 <- grepl("\\bnoOPIOD_No_NAS-Contr\\b", colnames(myNorm.PBC))
Group3 <- myNorm.PBC[,Group3]




## Folder 1  ## Group1 Vs Group2

folder1 <- as.data.frame(read_excel("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/I/Folder_1_7.xlsx"))


sample_set <- merge(Group1,Group2, by="row.names")
rownames(sample_set) <- sample_set$Row.names
sample_set <- sample_set[,c(-1)]

## select probes 899/982 probes 
sample_set <- sample_set[rownames(sample_set) %in% folder1$TargetID[grepl("cg",folder1$TargetID)],] 





## PCA Plot
pca <- prcomp(t(sample_set), scale. = TRUE)
pd <- ifelse(grepl("OPIOD_NAS-Aff", colnames(sample_set)),"OPIOD_NAS-Aff","OPIOD_No_NAS-Contr")
color <- as.factor(pd)
PCi<-data.frame(pca$x,Sample=color)

## 2D plot
ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("OPIOD_NAS-Aff"="red4","OPIOD_No_NAS-Contr"="blue4"))+ #your colors here
  theme_classic()

ggsave("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/I/2D_PCA_selected.pdf", dpi = 600, width = 12, height = 8)




################### 3D PCA plot ######################
groups <- levels(color)
#group.color <- c("noOPIOD"="blue4","OPIOD_NAS"="red4","OPIOD_No_NAS"="green4")
#colors <- c(rep("red4", 23), rep("blue4", 24))
#pca$pcolor <- colors


pca$pcolor <- pd

pca$pcolor <- ifelse(grepl("OPIOD_NAS-Aff", pca$pcolor), "red4", "blue4")

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 60)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("bottomleft", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("red4","blue4"), pch = 19)
dev.print(pdf, 'C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/I/3D_PCA_withXY.pdf', width = 12, height = 12)


########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80)
legend("bottomleft", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("red4","blue4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80),labels = rownames(pca$x), cex= 0.6, col = pca$pcolor)
dev.print(pdf, 'C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/I/3D_PCA_withXY_withName.pdf', width = 12, height = 12)




#####

rm(list=ls()[! ls() %in% c('myNorm.PBC','Group1','Group2','Group3')])


## Folder 2  ## (Group1 + Group2)  Vs Group3

folder2 <- as.data.frame(read_excel("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/II/Folder_2_8.xlsx"))

# not necessary but helps reorder
#sample_set <- Reduce(function(x, y) merge(x, y, by ="row.names"), list(Group1, Group2, Group3))
#sample_set <- reshape::merge_all(list(Group1, Group2, Group3), by= "row.names")

sample_set <- merge(Group1,Group2, by="row.names")
rownames(sample_set) <- sample_set$Row.names
sample_set <- sample_set[,c(-1)]
sample_set <- merge(sample_set, Group3, by ="row.names")
rownames(sample_set) <- sample_set$Row.names
sample_set <- sample_set[,c(-1)]

##  
sample_set <- sample_set[rownames(sample_set) %in% folder2$TargetID[grepl("cg",folder2$TargetID)],] 





## PCA Plot
pca <- prcomp(t(sample_set), scale. = TRUE)
pd <- ifelse(grepl("\\bOPIOD_NAS-Aff\\b|\\bOPIOD_No_NAS-Contr\\b", colnames(sample_set)),"OPIOD Used","No OPIOD Used")
color <- as.factor(pd)
PCi<-data.frame(pca$x,Sample=color)

## 2D plot
ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("OPIOD Used" = "red4",
                                "No OPIOD Used" ="blue4"))+ #your colors here
  theme_classic()

ggsave("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/II/2D_PCA_selected.pdf", dpi = 600, width = 12, height = 8)




################### 3D PCA plot ######################
groups <- levels(color)
#group.color <- c("noOPIOD"="blue4","OPIOD_NAS"="red4","OPIOD_No_NAS"="green4")
#colors <- c(rep("red4", 23), rep("blue4", 24))
#pca$pcolor <- colors


pca$pcolor <- pd

pca$pcolor <- ifelse(grepl("^\\bOPIOD\\sUsed\\b$", pca$pcolor), "red4", "blue4")

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 60)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("bottomleft", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
dev.print(pdf, 'C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/II/3D_PCA_withXY.pdf', width = 12, height = 12)


########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80)
legend("bottomleft", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80),labels = rownames(pca$x), cex= 0.6, col = pca$pcolor)
dev.print(pdf, 'C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/II/3D_PCA_withXY_withName.pdf', width = 12, height = 12)





## Folder 3

#####

rm(list=ls()[! ls() %in% c('myNorm.PBC','Group1','Group2','Group3')])


folder3 <- as.data.frame(read_excel("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/III/Folder_3_8.xlsx"))



sample_set <- merge(Group1,Group3, by="row.names")
rownames(sample_set) <- sample_set$Row.names
sample_set <- sample_set[,c(-1)]

## select probes 899/982 probes 
sample_set <- sample_set[rownames(sample_set) %in% folder3$TargetID[grepl("cg",folder3$TargetID)],] 





## PCA Plot
pca <- prcomp(t(sample_set), scale. = TRUE)
pd <- ifelse(grepl("OPIOD_NAS-Aff", colnames(sample_set)),"OPIOD_NAS-Aff","noOPIOD_No_NAS-Contr")
color <- as.factor(pd)
PCi<-data.frame(pca$x,Sample=color)

## 2D plot
ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("OPIOD_NAS-Aff"="red4","noOPIOD_No_NAS-Contr"="blue4"))+ #your colors here
  theme_classic()

ggsave("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/III/2D_PCA_selected.pdf", dpi = 600, width = 12, height = 8)




################### 3D PCA plot ######################
groups <- levels(color)
#group.color <- c("noOPIOD"="blue4","OPIOD_NAS"="red4","OPIOD_No_NAS"="green4")
#colors <- c(rep("red4", 23), rep("blue4", 24))
#pca$pcolor <- colors


pca$pcolor <- pd

pca$pcolor <- ifelse(grepl("OPIOD_NAS-Aff", pca$pcolor), "red4", "blue4")

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 60)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("bottomleft", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
dev.print(pdf, 'C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/III/3D_PCA_withXY.pdf', width = 12, height = 12)


########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80)
legend("bottomleft", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80),labels = rownames(pca$x), cex= 0.6, col = pca$pcolor)
dev.print(pdf, 'C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/III/3D_PCA_withXY_withName.pdf', width = 12, height = 12)




#####

#Folder 2: Probes taken from DiffMeth_AUC_Group1_and_2_to_Group3.csv



rm(list=ls()[! ls() %in% c('myNorm.PBC','Group1','Group2','Group3')])


## Folder 2_1  ## (Group1 + Group2)  Vs Group3

folder2 <- read.csv("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/DiffMeth_AUC_Group1_and_2_to_Group3.csv", header = TRUE, stringsAsFactors = FALSE)


# not necessary but helps reorder
#sample_set <- Reduce(function(x, y) merge(x, y, by ="row.names"), list(Group1, Group2, Group3))
#sample_set <- reshape::merge_all(list(Group1, Group2, Group3), by= "row.names")

sample_set <- merge(Group1,Group2, by="row.names")
rownames(sample_set) <- sample_set$Row.names
sample_set <- sample_set[,c(-1)]
sample_set <- merge(sample_set, Group3, by ="row.names")
rownames(sample_set) <- sample_set$Row.names
sample_set <- sample_set[,c(-1)]

##  
sample_set <- sample_set[rownames(sample_set) %in% folder2$X[grepl("cg",folder2$X)],] 





## PCA Plot
pca <- prcomp(t(sample_set), scale. = TRUE)
pd <- ifelse(grepl("\\bOPIOD_NAS-Aff\\b|\\bOPIOD_No_NAS-Contr\\b", colnames(sample_set)),"OPIOD Used","No OPIOD Used")
color <- as.factor(pd)
PCi<-data.frame(pca$x,Sample=color)

## 2D plot
ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("OPIOD Used" = "red4",
                                "No OPIOD Used" ="blue4"))+ #your colors here
  theme_classic()

ggsave("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/II/2D_PCA_selected_DiffMeth.pdf", dpi = 600, width = 12, height = 8)




################### 3D PCA plot ######################
groups <- levels(color)
#group.color <- c("noOPIOD"="blue4","OPIOD_NAS"="red4","OPIOD_No_NAS"="green4")
#colors <- c(rep("red4", 23), rep("blue4", 24))
#pca$pcolor <- colors


pca$pcolor <- pd

pca$pcolor <- ifelse(grepl("^\\bOPIOD\\sUsed\\b$", pca$pcolor), "red4", "blue4")

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 60)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("bottomleft", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
dev.print(pdf, 'C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/II/3D_PCA_DiffMeth.pdf', width = 12, height = 12)


########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80)
legend("bottomleft", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80),labels = rownames(pca$x), cex= 0.6, col = pca$pcolor)
dev.print(pdf, 'C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/II/3D_PCA_DiffMeth_withName.pdf', width = 12, height = 12)


####


rm(list=ls()[! ls() %in% c('myNorm.PBC','Group1','Group2','Group3')])


## Folder 2_1  ## (Group1 + Group2)  Vs Group3

folder2 <- read.csv("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/DiffMeth_AUC_Group1_and_2_to_Group3.csv", header = TRUE, stringsAsFactors = FALSE)


folder2 <- folder2[order(abs(folder2$deltaBeta), decreasing = TRUE),]

folder2 <- folder2[1:100,]


# not necessary but helps reorder
#sample_set <- Reduce(function(x, y) merge(x, y, by ="row.names"), list(Group1, Group2, Group3))
#sample_set <- reshape::merge_all(list(Group1, Group2, Group3), by= "row.names")

sample_set <- merge(Group1,Group2, by="row.names")
rownames(sample_set) <- sample_set$Row.names
sample_set <- sample_set[,c(-1)]
sample_set <- merge(sample_set, Group3, by ="row.names")
rownames(sample_set) <- sample_set$Row.names
sample_set <- sample_set[,c(-1)]

##  
sample_set <- sample_set[rownames(sample_set) %in% folder2$X[grepl("cg",folder2$X)],] 





## PCA Plot
pca <- prcomp(t(sample_set), scale. = TRUE)
pd <- ifelse(grepl("\\bOPIOD_NAS-Aff\\b|\\bOPIOD_No_NAS-Contr\\b", colnames(sample_set)),"OPIOD Used","No OPIOD Used")
color <- as.factor(pd)
PCi<-data.frame(pca$x,Sample=color)

## 2D plot
ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("OPIOD Used" = "red4",
                                "No OPIOD Used" ="blue4"))+ #your colors here
  theme_classic()

ggsave("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/II/2D_PCA_100_DiffMeth.pdf", dpi = 600, width = 12, height = 8)




################### 3D PCA plot ######################
groups <- levels(color)
#group.color <- c("noOPIOD"="blue4","OPIOD_NAS"="red4","OPIOD_No_NAS"="green4")
#colors <- c(rep("red4", 23), rep("blue4", 24))
#pca$pcolor <- colors


pca$pcolor <- pd

pca$pcolor <- ifelse(grepl("^\\bOPIOD\\sUsed\\b$", pca$pcolor), "red4", "blue4")

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 60)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("bottomleft", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
dev.print(pdf, 'C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/II/3D_PCA_100_DiffMeth.pdf', width = 12, height = 12)


########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80)
legend("bottomleft", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80),labels = rownames(pca$x), cex= 0.6, col = pca$pcolor)
dev.print(pdf, 'C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/PCA_Plot/II/3D_PCA_100_DiffMeth_withName.pdf', width = 12, height = 12)
