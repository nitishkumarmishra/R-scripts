Meth.TSS1500 <- TSS1500.order.aggregate[ TSS1500.order.aggregate$Group.1 %in% rownames(Expr),]
Expr.TSS1500 <- Expr[rownames(Expr)%in%Meth.TSS1500$Group.1,]
rownames(Meth.TSS1500) <- Meth.TSS1500$Group.1
Meth.TSS1500$Group.1 <- NULL
dim(Meth.TSS1500)
dim(Expr.TSS1500)
Expr.TSS1500 <- as.matrix(Expr.TSS1500)
Meth.TSS1500 <- as.matrix(Meth.TSS1500)
#Expr.TSS1500 <- log(Expr.TSS1500)
ll <- mapply(function(x,y)cor.test(Meth.TSS1500[x,],Expr.TSS1500[y,], method = "pearson", alternative = "greater"),
             1:nrow(Meth.TSS1500),
             1:nrow(Meth.TSS1500),
             SIMPLIFY=FALSE)
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
TSS1500.p <- cbind(p.value, cor.value)
rownames(TSS1500.p) <- rownames(Meth.TSS1500)
TSS1500.p <- data.frame(TSS1500.p)
adjust <- mt.rawp2adjp(TSS1500.p$p.value, proc=(c("Bonferroni")))
adjust.order <- adjust$adjp[order(adjust$index),]
adjust.order.1 <- data.frame(adjust.order)
final.TSS1500 <- cbind(TSS1500.p, adjust.order.1)
TSS1500.adjP.0.05 <- final.TSS1500[(final.TSS1500$rawp <= 0.05) & (final.TSS1500$Bonferroni <= 0.05),]
dim(TSS1500.adjP.0.05)
