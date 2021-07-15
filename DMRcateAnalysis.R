
library(pacman)
### p_load command can load several package at a time
p_load(ELMER, TCGAbiolinks, impute, matrixStats, ChAMP, limma, multtest, DESeq2, edgeR, DMRcate, missMethyl)
#p_loaded() ### To check the loaded packages
#p_unload(negate = TRUE)## Unload all packages
############# load TCGA level3 data ############
getTCGA("CHOL", Meth = TRUE, RNA = TRUE, Clinic = TRUE,   basedir="~", RNAtype = "gene", Methfilter = 0.2)
load("CHOL_clinic.rda")
load("CHOL_meth.rda")
load("CHOL_RNA.rda")
####################################
colnames(Meth) <- substr(colnames(Meth),1,16)
numNAs <- rowSums(is.na(Meth))
Meth <- Meth[!(numNAs > dim(Meth)[2]*0.25),]
Meth.impue <- impute.knn(Meth, k = 15, rowmax = 0.2, colmax = 0.2, rng.seed = 12345)
Meth.impue <- Meth.impue$data
####################################
myNorm <- champ.norm(beta = Meth.impue, fromIDAT = F, methValue = "B", norm = "BMIQ", fromFile = FALSE, betaFile = Meth.impue, filterXY = FALSE, plotBMIQ = TRUE, resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
Meth.impue.BMIQ <- myNorm$beta
Tumor.BMIQ <- Meth.impue.BMIQ[,which(substr(colnames(Meth.impue.BMIQ),14,15)=="01")]
Normal.BMIQ <- Meth.impue.BMIQ[,which(substr(colnames(Meth.impue.BMIQ),14,15)=="11")]
BMIQ.Meth <- cbind(Tumor.BMIQ, Normal.BMIQ)
#BMIQ.Meth.NoSNP <- rmSNPandCH(BMIQ.Meth, dist = 2, mafcut = 0.05)
BMIQ.Meth <- BMIQ.Meth[-grep("rs", rownames(BMIQ.Meth)),]
BMIQ.Meth.M <- logit2(BMIQ.Meth)
c1 <- ncol(Tumor.BMIQ)
c2 <- ncol(Normal.BMIQ)
###########DMRcate ##################
BMIQ.Meth.M.NoSNP <- rmSNPandCH(BMIQ.Meth.M, dist = 2, mafcut = 0.05)
#BMIQ.Meth.M.NoSNP.1 <- BMIQ.Meth.M.NoSNP[-grep("rs", rownames(BMIQ.Meth.M.NoSNP)),]
#### Design contrast matrix #########
#design = cbind(Tumor = c(rep(1,  c1), rep(0, c2)), Normal = c(rep(0, c1), rep(1, c2)))
#contrast.matrix = makeContrasts(TumorVsNormal = Tumor-Normal, levels = design)
groups <- as.factor(c(rep("Tumor",c1),rep("Normal",c2)))
design<-model.matrix(~0+groups)
colnames(design)=levels(groups)
contrast.matrix <- makeContrasts(Tumor-Normal, levels = design)
######### DMRcate anaysis ###########
myannotation <- cpg.annotate("array",BMIQ.Meth.M.NoSNP, analysis.type = "differential", design = design, contrasts = TRUE, cont.matrix = contrast.matrix, coef = "Tumor - Normal", fdr = 0.01)
#myannotation <- cpg.annotate("array",BMIQ.Meth.M.NoSNP, analysis.type = "differential", design = design, contrasts = TRUE, cont.matrix = contrast.matrix, coef = "Tumor - Normal", fdr = 0.01)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, p.adjust.method = "BH")
data(dmrcatedata)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
groups <- c(Tumor="magenta", Normal="forestgreen")
type <- as.factor(c(rep("Tumor", c1), rep("Normal", c2)))
cols <- groups[as.character(type)]
samps <- c(1:36, 36+(1:9))
#DMR.plot(ranges=results.ranges, dmr=1, CpGs=BMIQ.Meth, phen.col=cols, genome="hg19", samps=samps)
DMR.plot(ranges=results.ranges, dmr=1, CpGs=BMIQ.Meth, phen.col=cols, genome="hg19", samps=samps, showSampleNames = TRUE, cex.sampleNames = 0.8, separator = 1)
###################################################

###differential variability (DiffVar) from missMethyl ######
fitvar.contr <- varFit(BMIQ.Meth.M, design=design, coef=c(1,2))
fitvar.contr <- contrasts.varFit(fitvar.contr,contrasts=contrast.matrix)
summary(decideTests(fitvar.contr))
#topVar(fitvar.contr, number = 100)
VarDiff <- topVar(fitvar.contr,coef=1, sort = TRUE, number = 1000)
write.csv(VarDiff, file = "VarDiffmissMethyl.csv")

############## Conver Grange file in matrix #######
df1 <- data.frame(seqnames=seqnames(results.ranges),
                  starts=start(results.ranges)-1,
                  ends=end(results.ranges),
                  names=c(rep(".", length(results.ranges))),
                  strand=strand(results.ranges)
)
df <- mcols(results.ranges)
df <- mcols(results.ranges)