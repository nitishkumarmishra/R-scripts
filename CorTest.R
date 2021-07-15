# N <- 3751 * 1900
# cer.m <- matrix(1:N,ncol=1900)
# par.m <- matrix(1:N+rnorm(N),ncol=1900)
#aggregate(TSS200.order[,4:181], by=list(TSS200.order$GeneSymbol), FUN=median)
Meth.TSS1st <- TSS1stExon.order.aggregate[ TSS1stExon.order.aggregate$Group.1%in%rownames(Expr),]
rownames(Meth.TSS1st) <- Meth.TSS1st$Group.1
Meth.TSS1st$Group.1 <- NULL
Expr.TSS1st <- Expr[rownames(Expr)%in%Meth.TSS1st$Group.1,]
Meth.TSS1st <- as.matrix(Meth.TSS1st)
Expr.TSS1st <- as.matrix(Expr.TSS1st)

ll <- mapply(function(x,y)cor.test(Meth.TSS1st[x,],Expr.TSS1st[y,], method = "pearson", alternative = "greater"),
             1:nrow(Meth.TSS1st),
             1:nrow(Meth.TSS1st),
             SIMPLIFY=FALSE)
#########Same command with split #############
# ll <- mapply(cor.test,
#              split(par.m,row(par.m)),
#              split(cer.m,row(cer.m)),
#              SIMPLIFY=FALSE)
##############################################
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
TSS1st.p <- cbind(p.value, cor.value)
rownames(TSS1st.p) <- rownames(Meth.TSS1st)
TSS1st.p <- data.frame(TSS1st.p)
adjust <- mt.rawp2adjp(TSS1st.p$p.value, proc=(c("Bonferroni")))
adjust.order <- adjust$adjp[order(adjust$index),]
adjust.order.1 <- data.frame(adjust.order)
final.TSS1st <- cbind(TSS1st.p, adjust.order.1)
TSS1st.adjP.0.5 <- final.TSS1st[(final.TSS1st$rawp <= 0.05) & (final.TSS1st$Bonferroni <= 0.05),]