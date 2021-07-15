#### R program for pathway analysis by using clusterProfiler
detach("package:dplyr", unload=TRUE)
gene <- diffExp$geneID
geneList <- diffExp$log2FoldChange
name(geneList) <- diffExp$geneID
geneList <- sort(geneList, decreasing = TRUE)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 5,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
summary(kk2)
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa00982",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

gmtfile <- system.file("extdata", "c5.cc.v5.1.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
egmt <- enricher(gene, TERM2GENE=c5)
summary(egmt)


gmtfile <- system.file("extdata", "c6.all.v5.1.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
egmt <- enricher(gene, TERM2GENE=c5)
summary(egmt)
