# R code for Dr. Hamid Band
## In this code I am using cBioportal dta, which is normalized data from GDAC Firehose
BRCA_GeneExp <- read.csv("data_RNA_Seq_v2_expression_median.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
BRCA_GeneExp$Entrez_Gene_Id <- NULL
save(BRCA_GeneExp, file = "BRCA_cBioPortal_exp.Rda")
##########################################################
BRCA_GeneExp <- round(data.matrix(BRCA_GeneExp),0)
BRCA_GeneExp.log <- log(BRCA_GeneExp+1)
Cancer_BRCA_GeneExp.log <- BRCA_GeneExp.log[,which(substr(colnames(BRCA_GeneExp), 14,15)=="01")]

############## File for genes ##########
CTSB <- Cancer_BRCA_GeneExp.log[which(rownames(Cancer_BRCA_GeneExp.log)=="CTSB"),]
CTSL <- Cancer_BRCA_GeneExp.log[which(rownames(Cancer_BRCA_GeneExp.log)=="CTSL"),]
MZF1 <- Cancer_BRCA_GeneExp.log[which(rownames(Cancer_BRCA_GeneExp.log)=="MZF1"),]
######## Spearman correlation ##########
cor.test(MZF1, CTSB, method = "spearman")
cor.test(MZF1, CTSL, method = "spearman")
cor.test(CTSB, CTSL, method = "spearman")
######## Perason correlation ###########
cor.test(MZF1, CTSB, method = "pearson")
cor.test(MZF1, CTSL, method = "pearson")
cor.test(CTSB, CTSL, method = "pearson")
########### Plot #######################
######### In these plots MZF1 is on X-axis ##################
pdf("AllPlot.pdf", width = 10, height = 10)
par(mfrow=c(1,3))
plot(MZF1, CTSB, pch = 16, cex = 1, col = "black", main = "MZF1 Vs CTSB")
abline(lm(CTSB~MZF1))

plot(MZF1, CTSL, pch = 16, cex = 1, col = "black", main = "MZF1 Vs CTSL")
abline(lm(CTSL~MZF1))

plot(CTSB, CTSL, pch = 16, cex = 1, col = "black", main = "CTSB Vs CTSL")
abline(lm(CTSL~CTSB))
dev.off()
#######################################

save.image(file = "cBioPortla.RData")
######### In these plots MZF1 is on X-axis ##################
pdf("MZF1_vs_CTSB.pdf", width = 7, height = 7)
r <- cor.test(MZF1, CTSB,method = "spearman")$estimate[[1]]
my.p <- cor.test(MZF1, CTSB,method = "spearman")$p.value ## Fr
plot(MZF1, CTSB, pch = 16, cex = 1, col = "black", abline(lm(CTSB~MZF1), lwd=2))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R) == MYVALUE),
                   list(MYVALUE = format(r,dig=3)))[2]
rp[2] = substitute(expression(italic(P) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, cex = 0.8)
dev.off()

pdf("MZF1_vs_CTSL.pdf", width = 7, height = 7)
r <- cor.test(MZF1, CTSL,method = "spearman")$estimate[[1]]
my.p <- cor.test(MZF1, CTSL,method = "spearman")$p.value ## Fr
plot(MZF1, CTSL, pch = 16, cex = 1, col = "black", abline(lm(CTSL~MZF1), lwd=2))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R) == MYVALUE),
                   list(MYVALUE = format(r,dig=3)))[2]
rp[2] = substitute(expression(italic(P) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, cex = 0.8)
dev.off()

pdf("CTSB_vs_CTSL.pdf", width = 7, height = 7)
r <- cor.test(CTSB, CTSL,method = "spearman")$estimate[[1]]
my.p <- cor.test(CTSB, CTSL,method = "spearman")$p.value ## my.p <- 2e-16 
plot(CTSL, CTSB, pch = 16, cex = 1, col = "black", abline(lm(CTSB~CTSL), lwd=2))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R) == MYVALUE),
                   list(MYVALUE = format(r,dig=3)))[2]
rp[2] = substitute(expression(italic(P) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, cex = 0.8)
dev.off()

#################### Pearsn's correlation plot ################
pdf("MZF1_vs_CTSB_Pearson.pdf", width = 7, height = 7)
r <- cor.test(MZF1, CTSB,method = "pearson")$estimate[[1]]
my.p <- cor.test(MZF1, CTSB,method = "pearson")$p.value ## Fr
plot(MZF1, CTSB, pch = 16, cex = 1, col = "black", abline(lm(CTSB~MZF1), lwd=2))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R) == MYVALUE),
                   list(MYVALUE = format(r,dig=3)))[2]
rp[2] = substitute(expression(italic(P) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, cex = 0.8)
dev.off()

pdf("MZF1_vs_CTSL_Pearson.pdf", width = 7, height = 7)
r <- cor.test(MZF1, CTSL,method = "pearson")$estimate[[1]]
my.p <- cor.test(MZF1, CTSL,method = "pearson")$p.value ## Fr
plot(MZF1, CTSL, pch = 16, cex = 1, col = "black", abline(lm(CTSL~MZF1), lwd=2))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R) == MYVALUE),
                   list(MYVALUE = format(r,dig=3)))[2]
rp[2] = substitute(expression(italic(P) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, cex = 0.8)
dev.off()

pdf("CTSB_vs_CTSL_Pearson.pdf", width = 7, height = 7)
r <- cor.test(CTSB, CTSL,method = "pearson")$estimate[[1]]
my.p <- cor.test(CTSB, CTSL,method = "pearson")$p.value ## my.p <- 2e-16 
plot(CTSL, CTSB, pch = 16, cex = 1, col = "black", abline(lm(CTSB~CTSL), lwd=2))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R) == MYVALUE),
                   list(MYVALUE = format(r,dig=3)))[2]
rp[2] = substitute(expression(italic(P) == MYOTHERVALUE),
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, cex = 0.8)
dev.off()