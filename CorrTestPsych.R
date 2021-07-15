library("psych")
mat<-matrix(rnorm(60), nrow=12)
row.names(mat)<-rep(letters[1:c(nrow(mat)/3)], each=3)
#aggregate(TSS200.order[,4:181], by=list(TSS200.order$GeneSymbol), FUN=median)
aggregate(mat, by=list(row.names(mat)), FUN=median)
mat<-matrix(rnorm(60), nrow=12)
mat1<-matrix(rnorm(60), nrow=12)
colnames(mat) <- c("A", "B", "C", "D", "E")
colnames(mat1) <- c("A", "B", "C", "D", "E")                  
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat, diag = TRUE)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
res <- corr.test(mat,mat1,adjust="bonferroni",alpha=.05)
tmp <- flattenCorrMatrix(res$r, res$p)
tmp1 <- tmp[tmp$row==tmp$column,]
