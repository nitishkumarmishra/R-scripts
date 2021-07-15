k<-length(y)
for(i in 1:k) 
{ 
  nam <- paste("A", i, sep = "")
  assign(nam, listprobes[[i]])
}