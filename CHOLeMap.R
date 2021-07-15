### I have to run this tool in /storage/gudalab/nmishra/ICGC-data/TCGA-PAAD-Data/RNASeqV2/eMap/CHOL-eMAP
library(eMap) # eQTL analysis tool
library(quantsmooth)### Cytoband plot
res <- read.csv("refseq_genes_hg19_finished.txt", header=TRUE, sep="\t")### Refseq file with tss
probeinfo <- read.csv("Probeinfo.txt", header=TRUE, sep="\t")# Illumina probe annotation
load("GeneExp.Rda")### CHOL TCGA level3 RSEM expected count dev
load("BMIQ_Meth.Rda") #### BMIQ normalized CHOL level3 beta value data
expr <- round(geneExp)
expr <- log2(expr+1)
expr$gene <- sapply(strsplit(rownames(expr),"\\|"),'[[',1)
exprMerge <- merge(expr, res, by.x="gene", by.y="name2")
exprMerge$chrom <- gsub("chr","", exprMerge$chrom)

exprMerge$chrom <- gsub("chr","", exprMerge$chrom)
expr<-expr[order(expr$Chr,expr$TSS),]
me<-as.matrix(expr[,-(1:4)]);rownames(me)<-expr[,1]
eChr<-as.integer(expr$Chr)
ePos<-as.numeric(expr$TSS)
meth<-cbind(IlmnID=rownames(beta),data.frame(beta)) # Where "beta" is a matrix of methylation values.
meth<-merge(probeinfo,meth,by.x="Probe",by.y="IlmnID",all.x=F,all.y=F)
meth<-meth[order(meth$Chr,meth$Pos),]
mm<-as.matrix(meth[,-(1:7)]);rownames(mm)<-meth[,1]
mm<-mm[,colnames(mm)%in%colnames(me)]
me<-me[,colnames(me)%in%colnames(mm)]
me<-me[,order(colnames(me))]
mm<-mm[,order(colnames(mm))]
mChr<-meth$Chr
 mChr[mChr=="X"]<-23 #mChr <- gsub("X", "23", mChr, perl = TRUE) #I used gsub
 mChr[mChr=="Y"]<-24 #mChr <- gsub("X", "23", mChr, perl = TRUE)
 mChr<-as.integer(mChr)
mPos<-as.numeric(meth$Pos)
eMap1(me=me,mm1=mm,output.tag="output",p.cut=1,cis.only=T,cis.distance=1e5,eChr=eChr,ePos=ePos,mChr=mChr,mPos=mPos)
eqtl<-read.table("output_eqtl.txt",header=T,sep="\t",na.strings="NA",row.names=NULL)
eqtl<-eqtl[,colnames(eqtl)%in%c("Gene_ID","Marker_ID","b1","b1_p")]
eqtl$bonf_p<-p.adjust(eqtl$b1_p,method="bonferroni")
eannot<-cbind(data.frame(nr=1:nrow(me)),Gene=rownames(me),eChr=eChr,ePos=ePos)
mannot<-cbind(data.frame(nr=1:nrow(mm)),Probe=rownames(mm),mChr=mChr,mPos=mPos)
eqtl<-merge(eqtl,eannot,all.x=T,all.y=F,by.x="Gene_ID",by.y="nr")
eqtl<-merge(eqtl,mannot,all.x=T,all.y=F,by.x="Marker_ID",by.y="nr")
#eqtl<-merge(eqtl,expr[,colnames(expr)%in%c("Gene","Strand")],all.x=T,all.y=F,by.x="Gene",by.y="Gene")
eqtl<-merge(eqtl,me2[,colnames(me2)%in%c("Gene","Strand")],all.x=T,all.y=F,by.x="Gene",by.y="Gene")
eqtl$dist<-(eqtl$mPos-eqtl$ePos)*eqtl$Strand
saveRDS(eqtl,"output_eqtl_annot.rds")
eqtl<-eqtl[eqtl$bonf_p<0.05,]
saveRDS(eqtl,"output_eqtl_annot_bonf.rds")

pick<-eqtl$b1>0
logp<-(-log10(as.numeric(eqtl$b1_p[pick])))
mpos<-as.numeric(eqtl$mPos[pick])
tss<-as.numeric(eqtl$ePos[pick])
strand<-as.integer(eqtl$Strand[pick])
disttss<-(mpos-tss)*strand
pick<-disttss<1e5&disttss>-1e5
#pick<-disttss<1e4&disttss>-1e4
disttssf<-disttss[pick]
logpf<-logp[pick]
####### Negative correlations############
pick<-eqtl$b1<0
nlogp<-(-log10(as.numeric(eqtl$b1_p[pick])))
nmpos<-as.numeric(eqtl$mPos[pick])
ntss<-as.numeric(eqtl$ePos[pick])
nstrand<-as.integer(eqtl$Strand[pick])
ndisttss<-(nmpos-ntss)*nstrand
pick<-ndisttss<1e5&ndisttss>-1e5
ndisttssf<-ndisttss[pick]
nlogpf<-nlogp[pick]

pdf("Plot.pdf",width=10,height=8)
plot(ndisttssf,nlogpf,col="red",pch=20,cex=.6,
     ylim=c(min(c(logpf,nlogpf)),max(c(logpf,nlogpf))),
     xlab="Distance between CpG sites and TSS (bp)",
     ylab="-log(P-value)",
     xaxt="n"
)
#axis(1, xaxp=c(-100000, 100000, 5), las=2)
axis(1,at=seq(-1e5,1e5,length.out=5),labels=c("-100 000","-50 000", "0", "50 000","100 000"))
#axis(1,at=seq(-1e4,1e4,length.out=5),labels=c("-10 000","-5 000","0","5 000","10 000"))
legend("topleft", c("Red: Negative correlations","Blue: Positive correlations"), text.col=c("red","blue"),bty="n", col=c("red","blue"), lwd=1, lty=c(0,0),pch=c(20,20))
points(disttssf,logpf,col="blue",pch=20,cex=.6)
dev.off()
