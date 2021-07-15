a <- matrix(data=NA,nrow=194,ncol=2)
a[,1] <- colnames(beta.table)
region <- c(rep("Tumor", 184), rep("Normal", 10))
write.table(a, "TCGA-list.txt", row.names=FALSE, col.names=FALSE, sep = "\t")
beta.table <- t
### Replace . with _ inn colname of beta.table####
b <- colnames(beta.table)
#b1 <- gsub("\\.", "_", b)
b1 <- gsub("\\-", "_", b)
colnames(beta.table) <- b1
rm(b, b1)
######################################################
probes <- rownames(beta.table)
output.table <- data.frame(SiteID=probes, beta.table)
write.table(output.table, file="TCGA.txt", sep="\t", quote=F, row.names=F)
beta.file <- "TCGA.txt"
rm(probes)

library(COHCAP)
sample.file <- "TCGA-list.txt"
project.folder <- getwd()
project.name <- "TCGA-450k-DMR"

beta.table <- COHCAP.annotate(beta.file, project.name, project.folder, platform="450k-UCSC")
COHCAP.qc(sample.file, beta.table, project.name, project.folder)
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="Tumor", methyl.cutoff=0.5, unmethyl.cutoff = 0.1, delta.beta.cutoff = 0.2, pvalue.cutoff=0.01, fdr.cutoff=0.01, num.groups=2, paired=FALSE)
####COHCAP.avg.by.site does't work'
#filtered.avg.sites <- COHCAP.avg.by.site(filtered.sites, beta.table, project.name, project.folder, methyl.cutoff=0.5, unmethyl.cutoff = 0.1, delta.beta.cutoff = 0.2, pvalue.cutoff=0.01, fdr.cutoff=0.01, num.sites=3)
filtered.islands <- COHCAP.avg.by.island(sample.file, filtered.sites, beta.table, project.name, project.folder, ref="Tumor", methyl.cutoff=0.5, unmethyl.cutoff = 0.1, delta.beta.cutoff = 0.2, pvalue.cutoff=0.01, fdr.cutoff=0.01, num.groups=2, num.sites=3, plot.box=TRUE, paired=FALSE)
