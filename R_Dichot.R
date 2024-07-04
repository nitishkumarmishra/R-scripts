
metaData <- data.frame(error= numeric(0), Relate= integer(0), Dichot = integer(0), Y = integer(0))
nExp = 1000

error <- seq(5,250,2)
numTrial <- seq(1,nExp,1)

start <- proc.time()[3]
for (i in error ) {
    show(i)
   pDichot = list()
   pCont = list()
    for(exp in numTrial) {
      year <- round(runif(20, 1990,2030))
      words <-3.3*year
      words <- 6701 - words
      words <- words + rnorm(20, 0, i)
      split <- 1* (year <= 2010)
      dt <- matrix(c(year,words,split), ncol = 3)
      dt <- data.frame(dt)
      colnames(dt) <-c("year", "words", "split")
      modelC <- lm(words ~ year)
      modelD <- lm(words ~ split)
      pC <- summary(modelC)$coefficients[8]
      pD <- summary(modelD)$coefficients[8]
      pDichot <- c(pDichot, pD)
      pCont <- c(pCont, pC)

    };
  pDichot <-1*(pDichot <= 0.05)
  pCont <-1*(pCont <= 0.05)
  row <- c(i ,1,0,Reduce("+",pCont))
  metaData <- rbind(metaData,row)
  row <- c(i ,1,1,Reduce("+", pDichot))
  metaData <- rbind(metaData,row)

  
  pDichot = list()
  pCont = list()
  for(exp in numTrial) {
    year <- round(runif(20, 1990,2030))
    words <- rnorm(20, 100, i)
    split <- 1* (year <= 2010)
    dt <- matrix(c(year,words,split), ncol = 3)
    dt <- data.frame(dt)
    colnames(dt) <-c("year", "words", "split")
    modelC <- lm(words ~ year)
    modelD <- lm(words ~ split)
    pC <- summary(modelC)$coefficients[8]
    pD <- summary(modelD)$coefficients[8]
    pDichot <- c(pDichot, pD)
    pCont <- c(pCont, pC)
    
  };
  pDichot <-1*(pDichot <= 0.05)
  pCont <-1*(pCont <= 0.05)
  row <- c(i ,0,0,Reduce("+",pCont))
  metaData <- rbind(metaData,row)
  row <- c(i ,0,1,Reduce("+", pDichot))
  metaData <- rbind(metaData,row)

}
end <- proc.time()[3]
(end-start)/60


write.csv(metaData, '2.csv', row.names = FALSE)

