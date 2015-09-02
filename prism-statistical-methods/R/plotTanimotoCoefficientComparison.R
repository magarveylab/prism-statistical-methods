setwd("/Users/michaelskinnider/Desktop/ecfp6/v4/")

library(data.table)
library(ggplot2)
library(scales)

data <- read.csv("streptomyces.csv")
data$Software <- factor(data$Software, as.character(data$Software))
str(data)

plot <- ggplot(data=data, aes(x=Software, y=Tc, fill=Software)) + 
  geom_boxplot() + scale_y_continuous(limits = c(0, 1.0)) + theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank()) 
plot(plot)
