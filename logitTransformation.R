## R program for logit conversion of Beta value
logit <- function(p)
  log(p/(1-p))
logit2 <-   function(p)
  log2(p/(1-p))