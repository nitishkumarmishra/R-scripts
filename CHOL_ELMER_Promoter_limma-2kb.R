### R code for analysis of gene expression and methylation
library(ELMER)
load("CHOL_clinic.rda")
load("BMIQ_Meth.Rda")
load("CHOL_RNA.rda")
DiffMeth <- read.csv("DiffMethDelta.0.10.txt", header = TRUE, sep = ",", row.names = 1)
DiffMeth.0.25 <- DiffMeth[which(abs(DiffMeth$MeanDiff)>=0.25),]
colnames(GeneExp) <- substr(colnames(GeneExp), 1, 16)

save(GeneExp, file = "CHOL_RNA.rda")

probe <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, what="Locations")
probe <- GRanges(seqnames=probe$chr,
                 ranges=IRanges(probe$pos,width=1,names=rownames(probe)),
                 strand=probe$strand, name=rownames(probe))
#mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=TRUE, probeInfo=probe)
geneAnnot <- txs()
geneAnnot$GENEID <- paste0("ID",geneAnnot$GENEID)
geneInfo <- promoters(geneAnnot,upstream = 2000, downstream = 2000)
save(geneInfo,file="geneAnnot.rda")
mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=TRUE, probeInfo=probe, geneInfo=geneInfo)
getGeneInfo(mee)
Probe <- get.feature.probe(rm.chr=paste0("chr",c("X","Y", "M")), list(upstream = 2000, downstream = 2000), promoter = TRUE)
save(Probe,file="probeInfo_feature_distal.rda")


df1 <- data.frame(seqnames=names(Probe),
           starts=start(Probe),
           ends=end(Probe),
           names=c(rep(".", length(Probe))),
           strand=strand(Probe))
rownames(df1) <- df1$seqnames

mee <- fetch.mee(meth="BMIQ_Meth.Rda",exp="CHOL_RNA.rda", TCGA=TRUE,
                 probeInfo="probeInfo_feature_distal.rda",geneInfo="geneAnnot.rda")

########## Hypomethylated probes enhancer #####################
setwd("C:/Users/nitish.mishra/Desktop/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/ELMER/ELMER_Promoter_limma/Hypo")

#sig.diff <- get.diff.meth(mee, cores=detectCores()/2, diff.dir="hypo", pvalue = 0.01, sig.dif = 0.20)
#Sig.probes <- read.csv("getMethdiff.hypo.probes.significant.csv",stringsAsFactors=FALSE)[,1]
Sig.probes <- rownames(DiffMeth.0.25[which(DiffMeth.0.25$MeanDiff < 0.0),])
Sig.probes <- Sig.probes[Sig.probes%in%rownames(df1)]
nearGenes <-GetNearGenes(TRange=getProbeInfo(mee,probe=Sig.probes), geneAnnot=getGeneInfo(mee),geneNum = 2, cores=detectCores()/2)
Hypo.pair <-get.pair(mee=mee,probes=Sig.probes,nearGenes=nearGenes,
                     permu.dir="./permu",permu.size=200,Pe = 0.01
                     ,cores=detectCores()/2,label= "hypo")
Sig.probes.paired <- read.csv("getPair.hypo.pairs.significant.csv",stringsAsFactors=FALSE)[,1]
enriched.motif <-get.enriched.motif(probes=Sig.probes.paired, label="hypo", min.incidence = 10,lower.OR = 1.1)
load("getMotif.hypo.enriched.motifs.rda")
TF <- get.TFs(mee=mee, enriched.motif=enriched.motif, cores=detectCores()/2, label= "hypo")
#scatter.plot(mee,byProbe=list(probe=c("cg24874425"),geneNum=20),category="TN", save=TRUE)
#scatter.plot(mee,byPair=list(probe=c("cg24874425"),gene=c("ID51702")), category="TN", save=TRUE,lm_line=TRUE)

pair <- fetch.pair(pair="getPair.hypo.pairs.significant.withmotif.csv", probeInfo = "../probeInfo_feature_distal.rda", geneInfo = "../geneAnnot.rda")
#schematic.plot(pair=pair, byProbe="cg24874425",save=TRUE)
save.image(file = "ELMER_Hypo.RData")
################# Hyper methylated enhancer ##################
setwd("C:/Users/nitish.mishra/Desktop/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/ELMER/ELMER_Promoter_limma/Hyper")

#sig.diff <- get.diff.meth(mee, cores=detectCores()/2, diff.dir="hyper", pvalue = 0.01, sig.dif = 0.20)
#Sig.probes <- read.csv("getMethdiff.hyper.probes.significant.csv",stringsAsFactors=FALSE)[,1]
Sig.probes <- rownames(DiffMeth.0.25[which(DiffMeth.0.25$MeanDiff > 0.0),])
Sig.probes <- Sig.probes[Sig.probes%in%rownames(df1)]
nearGenes <-GetNearGenes(TRange=getProbeInfo(mee,probe=Sig.probes), geneAnnot=getGeneInfo(mee),geneNum = 2, cores=detectCores()/2)
Hyper.pair <-get.pair(mee=mee,probes=Sig.probes,nearGenes=nearGenes,
                      permu.dir="./permu",permu.size=200,Pe = 0.01,
                      cores=detectCores()/2,label= "hyper")
Sig.probes.paired <- read.csv("getPair.hyper.pairs.significant.csv",stringsAsFactors=FALSE)[,1]
enriched.motif <-get.enriched.motif(probes=Sig.probes.paired, label="hyper", min.incidence = 10,lower.OR = 1.1)
load("getMotif.hyper.enriched.motifs.rda")
TF <- get.TFs(mee=mee, enriched.motif=enriched.motif, cores=detectCores()/2, label= "hyper")
#scatter.plot(mee,byProbe=list(probe=c("cg24874425"),geneNum=20),category="TN", save=TRUE)
#scatter.plot(mee,byPair=list(probe=c("cg24874425"),gene=c("ID51702")), category="TN", save=TRUE,lm_line=TRUE)

pair <- fetch.pair(pair="getPair.hyper.pairs.significant.withmotif.csv", probeInfo = "../probeInfo_feature_distal.rda", geneInfo = "../geneAnnot.rda")
#schematic.plot(pair=pair, byProbe="cg24874425",save=TRUE)
save.image(file = "ELMER_Hyper.RData")
