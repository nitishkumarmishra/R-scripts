library(ggplot2)
library(reshape2)
library(scales)
stacked <- read.csv("Stacked1.txt", header = TRUE, row.names = 1, sep = "\t")
datm <- melt(cbind(stacked, GeneSubregion = rownames(stacked)), id.vars = c('GeneSubregion'))
ggplot(datm,aes(x = variable, y = value,fill = GeneSubregion)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_y_continuous(labels = percent_format())

