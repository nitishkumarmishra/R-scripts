library(psych)
#########Program to calculate correlation and p value for two matrix#######
####### Final out put give matrix #########################################
# row column cor p
# 1    A      A   0 1
# 3    B      B   0 1
# 6    C      C   0 1
# 10   D      D   0 1
# 15   E      E   0 1
##############################################################
mat<-matrix(rnorm(60), nrow=12)
mat1<-matrix(rnorm(60), nrow=12)
colnames(mat) <- c("A", "B", "C", "D", "E")
colnames(mat1) <- c("A", "B", "C", "D", "E")  
##############################################################
############## P.value for two matrix ########################
cor.test.p <- function(x,y){
  FUN <- function(x, y) cor.test(x, y, method = "pearson", alternative = "greater")[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(y), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}
##############################################################
############## Correlation for two matrix ####################
cor.test.cor <- function(x,y){
  FUN <- function(x, y) cor.test(x, y, method = "pearson", alternative = "greater")[["estimate"]]
  z <- outer(
    colnames(x), 
    colnames(y), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}
##############################################################
###########Function to generate matrix from P.value & correlation ################
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat, diag = TRUE)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
##############################################################
############## run for output ################################
tt <- cor.test.p(mat, mat1)
tt1 <- cor.test.cor(mat, mat1)
tmp <- flattenCorrMatrix(tt, tt1)
tmp1 <- tmp[tmp$row==tmp$column,]
