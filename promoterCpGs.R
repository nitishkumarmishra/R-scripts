library('org.Hs.eg.db')
library(Homo.sapiens)
library(biomaRt)
library(FDb.InfiniumMethylation.hg19)
#######################################################
genes <- read.csv("TSS-HGNC-name.csv") ### TSS-HGNC.CSV is the list of TCGA genes in RNASeqV2 level-3 data
yy<-as.character(genes$hgnc, na.rm=FALSE)
keytypes(Homo.sapiens)
x <- org.Hs.egSYMBOL2EG
mapped_genes <- mappedkeys(x)
gene <- intersect(yy, mapped_genes)
################################################
##List of HGNC symbol whose corresponding ENTREZ-ID available
txs <- transcriptsBy(Homo.sapiens, 'gene', col='GENEID')
hm450 <- get450k()
#use Entrez GeneID as a string, not numeric or factor
getProbes<-function(geneID){
  tmp <- org.Hs.egSYMBOL2EG[[geneID]]
  temp<-txs[[tmp]]
  if(is.null(temp)){
    probes<-NULL
  }else{
    upstream.probes<-names(subsetByOverlaps(hm450,flank(temp,1500)))
    downstream.probes<-names(subsetByOverlaps(hm450,flank(temp,-1500,start=TRUE)))
    probes<-unique(c(upstream.probes,downstream.probes))
  }
  return(probes)
}
gene1 <- gene[1:10]
sapply(gene1, getProbes)

############################ Final line to get list of probes for all genes ##########
##listprobe is character vector of all probes for yy, while lisprobes1 is list of CpGs for yy
listprobes <- unlist(sapply(yy, getProbes))
listprobes1 <-lapply(yy, getProbes)
