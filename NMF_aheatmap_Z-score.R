#mat <- exp.fusion.proteinCoding[common_names,]
#mat <- log2(mat+1)
load(file = "BRCA_Neetha_GSEA_Analysis.RData")
rna <- log2(BRCA.HTSeq.Counts[common_names,])
# get the index of the normal/control samples
Middle <- Neetha.freq.1.cBioPortal[Neetha.freq.1.cBioPortal$quantile=="Middle",]
#n_index <- which(substr(colnames(rna),14,14) == '1')
n_index <- which(substr(colnames(rna),1,12)%in%Middle$Row.names)
t_index <- which(substr(colnames(rna),14,14) == '0')

#z = [(value gene X in tumor Y)-(mean gene X in normal)]/(standard deviation X in normal)
rna_fusion <- rna[,colnames(exp.fusion.proteinCoding)]
mean_n <- rowMeans(rna[,n_index])
sd_n <- matrixStats::rowSds(as.matrix(rna[,n_index]))
z_rna <- rna_fusion - (mean_n)/ (sd_n)
z_rna <- na.omit(z_rna)
z_rna <- as.matrix(z_rna)
#z_rna[which(!is.finite(z_rna))] <- 0 ## replace all non-finite values with 0
z_rna <- z_rna[!rowSums(!is.finite(z_rna)),] ## remove all rows with non-finite values
z_rna1 <- z_rna[abs(rowMeans(z_rna)) >= 4,]

factors <- factor(c(rep("High", length(exp.High)), rep("Low", length(exp.Low))))
annoatation <- data.frame(RecurrentFusion=factors)

NMF::aheatmap(z_rna1, annCol = annoatation, color = "-RdBu:25")

# calculate z-scores "https://www.biostars.org/p/153013/"
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}
z_rna <- scal(rna_vm[,t_index],rna_vm[,n_index])


mat_scaled= t(apply(rna_fusion, 1, scale))
