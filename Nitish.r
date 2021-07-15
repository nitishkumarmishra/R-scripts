



# NB: There is no guarantee that this code will work for you precisely as it is now. But some tweeking on your side should make it work. 



# --- Illumina annotation
# http://support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html
data<-read.csv("HumanMethylation450_15017482_v1-2.csv",header=T,sep=",",na.strings="",row.names=NULL,quote="",skip=7)
data<-data[-(485578:nrow(data)),] #remove controls
probeinfo<-data[,colnames(data)%in%c("IlmnID","CHR","MAPINFO","Relation_to_UCSC_CpG_Island","Enhancer","Infinium_Design_Type","UCSC_CpG_Islands_Name")]
colnames(probeinfo)<-c("Probe","Design","Chr","Pos","CGI_name","Relation_CGI","Enhancer")
write.table(probeinfo, file="Probeinfo.txt", row.names=F, quote=F, sep="\t")




# --- Make table of genes containing gene name, TSS and strand. 

# Go to UCSC table browser (http://genome.ucsc.edu/cgi-bin/hgTables)
# assembly: hg19
# group: Genes and Gene Predictions
# track: RefSeq Genes
# table: refGenes
# output format: all fields from selected table
# output file: refseq_genes_hg19.txt

data<-read.table("refseq_genes_hg19.txt",header=T,sep="\t",na.strings="",row.names=NULL,quote="",comment.char="")
data<-data[,colnames(data)%in%c("name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","name2")]

# The following code removes alternative transcripts for each gene, so that each gene only has one TSS. 
res<-data[1,]
k<-0
for(g in unique(data$name2)){
	k<-k+1;print(k)
	d<-data[data$name2==g,]
	if(nrow(d)==1){res<-rbind(res,d)}
	if(nrow(d)>1){
		d<-d[grep("NM",d$name),]
		d$char<-as.numeric(lapply(as.character(d$name),nchar))
		d$num<-as.numeric(gsub("NM_","",d$name))
		d<-d[order(d$char,d$num),]
		d<-d[1,1:ncol(res)] # This approach will chose the NM number with 1) the fewest characters, and, if equal 2) the lowest number. This might not be ideal, but it't the best I could think of.
		res<-rbind(res,d)
	}
}
res<-res[-1,]
res<-res[!is.na(res$name2),]
write.table(res, file="refseq_genes_hg19_finished.txt", row.names=F, quote=F, sep="\t")

# Now you know the TSS of all genes, and from Illumina annotation you know the location of the probes


# --- eMap -- meth-expression analysis
library(eMap) # see references in paper for web page
expr<-read.table("...",header=T,sep="\t",na.strings="NA",row.names=NULL)
	# This file should contain 4 information columns with "Gene", "Chr", "TSS" and "Strand". The rest of the columns should be expression data of samples. 
	# or make the object in R.
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
 mChr[mChr=="X"]<-23
 mChr[mChr=="Y"]<-24
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
eqtl<-merge(eqtl,expr[,colnames(expr)%in%c("Gene","Strand")],all.x=T,all.y=F,by.x="Gene",by.y="Gene")
eqtl$dist<-(eqtl$mPos-eqtl$ePos)*eqtl$Strand
saveRDS(eqtl,"output_eqtl_annot.rds")
eqtl<-eqtl[eqtl$bonf_p<0.05,]
saveRDS(eqtl,"output_eqtl_annot_bonf.rds")

# Plot
# Positive correlations
pick<-eqtl$b1>0
logp<-(-log10(as.numeric(eqtl$b1_p[pick])))
mpos<-as.numeric(eqtl$mPos[pick])
tss<-as.numeric(eqtl$TSS[pick])
strand<-as.integer(eqtl$Strand[pick])
disttss<-(mpos-tss)*strand
pick<-disttss<1e5&disttss>-1e5
disttssf<-disttss[pick]
logpf<-logp[pick]
# Negative correlations
pick<-eqtl$b1<0
nlogp<-(-log10(as.numeric(eqtl$b1_p[pick])))
nmpos<-as.numeric(eqtl$mPos[pick])
ntss<-as.numeric(eqtl$TSS[pick])
nstrand<-as.integer(eqtl$Strand[pick])
ndisttss<-(nmpos-ntss)*nstrand
pick<-ndisttss<1e5&ndisttss>-1e5
ndisttssf<-ndisttss[pick] 
nlogpf<-nlogp[pick]

pdf("plot.pdf",width=10,height=8)
plot(ndisttssf,nlogpf,col="red",pch=20,cex=.6,
	ylim=c(min(c(logpf,nlogpf)),max(c(logpf,nlogpf))),
	xlab="Distance between methylation probe and TSS (bp)",
	ylab="-log p-value",
	xaxt="n"
	)
axis(1,at=seq(-1e5,1e5,length.out=5),labels=c("-100 000","-50 000","0","50 000","100 000"))
legend("topleft", c("Red: negative correlations","Blue: positive correlations"), text.col=c("red","blue"),bty="n")
points(disttssf,logpf,col="blue",pch=20,cex=.6)
dev.off()



# --- Quantsmooth: eQTL significance level
library(quantsmooth)
annot<-eqtl[,c(10,11)] # chromosome and position
colnames(annot)<-c("CHR","MapInfo")
annot$CHR[annot$CHR==23]<-"X"
annot$CHR[annot$CHR==24]<-"Y"
bmp("plot.bmp",height=800,width=800)
out<-prepareGenomePlot(annot,sexChromosomes=T,organism="hsa",units="bases",paintCytobands=F)
out<-data.frame(out)
logp<--log(eqtl$b1_p)
x<-out$MapInfo[eqtl$b1<0]
y<-(out$CHR+(logp/100))[eqtl$b1<0]
points(x=x,y=y,pch=20,cex=.4,col="red")
x<-out$MapInfo[eqtl$b1>0]
y<-(out$CHR+(logp/100))[eqtl$b1>0]
points(x=x,y=y,pch=20,cex=.5,col="blue")
dev.off()

















