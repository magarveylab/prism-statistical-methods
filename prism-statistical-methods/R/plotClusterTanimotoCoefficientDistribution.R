setwd("/Users/michaelskinnider/Desktop/")
data <- read.csv("meridamycin.csv")
library(ggplot2)
library(scales)
str(data)
plot <- ggplot(data=data, aes(x=Bin, y=Size)) + geom_bar(stat="identity") + theme_bw() +
  labs(x="Tanimoto coefficient",y="Predicted structures") + ggtitle("Meridamycin")
plot(plot)