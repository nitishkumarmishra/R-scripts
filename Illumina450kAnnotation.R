### This is R tool to make Illumina 450k probe annotation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
load("probe.features.Rda")## This is annotation data from chAMP
Illumina450k.probe.anno <- merge(ann450k, probe.features, by="row.names")
save(Illumina450k.probe.anno, file = "illIllumina450k.probe.anno.Rda")
