library(MethylMix)
data(METcancer)#### Cancer patients Methylation data
data(METnormal)### Normal paitents Methylation data
data(MAcancer) ### Gene expression data of cancer patients
MethylMixResults <- MethylMix(METcancer,METnormal,MAcancer, Parallel = TRUE)
#### This step fol plotting all genes #####
for (i in 1:length(MethylMixResults$MethylationDrivers)) 
{
  MethylMix_PlotModel(MethylMixResults$MethylationDrivers[i],METcancer,
  MethylMixResults,MAdata=0,METnormal)
}