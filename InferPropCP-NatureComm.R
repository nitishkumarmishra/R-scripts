#### InferPropCP.R
### Description: an R-script to infer the proportions of a priori known cell subtypes present in a sample representing a mixture of such cell-types. Inference proceeds via constrained projection (CP) and quadratic programming and thus requires the quadprog CRAN R package to be installed.
### Author: Andrew E Teschendorff (a.teschendorff@ucl.ac.uk)
### Date: 18th Sept. 2015
#### Copyright 2015 Andrew Teschendorff
#### Copyright permission: InferPropCP is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details (see <http://www.gnu.org/licenses/> )

#### LIBRARIES NEEDED
library(quadprog);

#### MAIN USER FUNCTION

#### INPUT ARGS:
### cent.m: a matrix of "centroids", i.e. representative molecular profiles, for a number of cell subtypes. rows label molecular features (e.g. genes, CpGs,...) and columns label the cell-type. IDs need to be provided as rownames and colnames, respectively. No missing values are allowed, and all values in this matrix should be positive or zero.
### data.m: a data matrix with rows labeling the molecular features (should use same ID as in cent.m) and columns labeling samples (e.g. primary tumour specimens). No missing values are allowed and all values should be positive or zero. Note that typically, nrow(data.m) >> nrow(cent.m), since the centroids are usually defined over a subset of discriminatory molecular features. So, these discriminatory features need to be identified beforehand, and centroids need to be constructed over these selected features.


### OUTPUT ARGS:
### westQP.m: a matrix representing the estimated cell proportions (weights) for each sample. Rows label samples, columns label the cell-types. So the entries in each row should add to 1.

### GUIDELINE for construction of cent.m: Given nCT cell types, the cent.m is best constructed using a supervised method like limma, comparing each cell-type in turn to all other cell-types, so nCT supervised analyses in total. Assuming that FDR thresholds allow identification of a fixed common number, ntop (typically in the range ~ 50-500), of cell-specific molecular features, then the cent.m should have a dimension of nCT*ntop X nCT (assuming no overlaps between the ntop discriminatory features of each cell-type). In practice, nrow(cent.m) < ntop*nCT. If FDR thresholds do not allow reliable identification of cell-type specific markers for every cell-type, consider relaxing FDR thresholds.

InferPropCP <- function(cent.m,data.m){
  
  require(quadprog);
  
  ### define D matrix
  nCT <- ncol(cent.m);
  D <- matrix(NA,nrow=nCT,ncol=nCT);
  for(j in 1:nCT){
    for(k in 1:nCT){
      D[j,k] <- 2*sum(cent.m[,j]*cent.m[,k])
    }
  }
  
  ### define constraints
  A.m <- matrix(0,nrow=nCT+1,ncol=nCT);
  A.m[1,] <- 1;
  for(i in 1:nCT){
    A.m[1+i,i] <- 1;
  }
  A.m <- t(A.m);
  b0.v <- c(1,rep(0,nCT));
  
  ### define d-vector and solve for each sample
  nS <- ncol(data.m);
  westQP.m <- matrix(NA,ncol=ncol(cent.m),nrow=nS);
  colnames(westQP.m) <- colnames(cent.m);
  rownames(westQP.m) <- colnames(data.m);
  
  match(rownames(cent.m),rownames(data.m)) -> map.idx;
  rep.idx <- which(is.na(map.idx)==FALSE);
  for(s in 1:nS){
    tmp.v <- data.m[,s];
    d.v <- as.vector(2*matrix(tmp.v[map.idx[rep.idx]],nrow=1) %*% cent.m[rep.idx,]);
    qp.o <- solve.QP(D,d.v,A.m,b0.v,meq=1);
    westQP.m[s,] <- qp.o$sol;
    print(s);
  }
  
  return(westQP.m);
  
}
