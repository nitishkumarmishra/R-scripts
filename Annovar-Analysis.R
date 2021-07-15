setwd("C:\\Users\\nitish.mishra\\Desktop\\GATK-Germline")
file<-list.files(full.names=TRUE)
k<-length(file)
for(i in 1:k) {
tmp <- read.csv("PACA-37N_final_snps.hg19_multianno.csv")
tmp2 <- tmp[which(tmp[,19]=="synonymous SNV"),]
write.csv(tmp2,"C:/Users/nitish.mishra/Desktop/outfile", sep = ",")
}
#"synonymous SNV"
#"nonsynonymous SNV"
#"frameshift insertion"
#"nonframeshift insertion"
#"frameshift deletion"
#"nonframeshift deletion")); 
#"stopgain SNV"
#"stoploss SNV"
#}
#rownames(sta)<-file
#colnames(sta) <- c("Number","synonymous SNV","nonsynonymous SNV","frameshift insertion","nonframeshift insertion","frameshift deletion","nonframeshift deletion","stopgain SNV","stoploss SNV","splicing","unknown")
#write.table(sta,"C:/Users/nitish.mishra/Desktop/MutationSummaryPCAll.csv", sep = ",")
