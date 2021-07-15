library(biomaRt)
# define biomart object
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# read in the file
genes <- read.csv("TSS-HGNC-name.csv")
# query biomart
results <- getBM(attributes = c("entrezgene", "hgnc_symbol"), filters = "hgnc_symbol", values = genes$hgnc, mart = mart)
