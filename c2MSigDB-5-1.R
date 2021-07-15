### R program to generate latest data from MSigDB5.1
### I can use this data as substitute of c2BroadSets
## Total 4746 gene compare to 3272 in c2BroadSets
library(GSEABase)
c2broad5.1.Symbol <- getGmt("C:/Users/nitish.mishra/Desktop/MSigDB/c2.all.v5.1.symbols.gmt", collectionType=BroadCollection(category="c2"), geneIdType=SymbolIdentifier())
c2broad5.1.Entrez <- getGmt("C:/Users/nitish.mishra/Desktop/MSigDB/c2.all.v5.1.entrez.gmt", collectionType=BroadCollection(category="c2"), geneIdType=SymbolIdentifier())
save.image(file = "c2MSigDB.5.1.RData")