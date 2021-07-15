########  Spearma's correlation by using permutation ##################
test <- function(x,y) suppressWarnings(cor.test(x, y, method="spearman")$p.value)
rho <- test(x,y)                                     # Test statistic
p <- replicate(10^3, test(sample(x, length(x)), sample(y, length(y))))   # Simulated permutation distribution

p.out <- sum(abs(p) < rho)    # Count of strict (absolute) exceedances
p.at <- sum(abs(p) == rho)    # Count of equalities, if any
(p.out + p.at /2) / length(p) # Proportion of exceedances: the p-value.
