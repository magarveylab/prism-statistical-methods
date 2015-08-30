# setwd("")

library(data.table)
library(ggplot2)
library(scales)

data <- read.csv("macrolide.csv")
data$Software <- factor(data$Software, as.character(data$Software))
str(data)

plot <- ggplot(data=data, aes(x=Software, y=Average, fill=Software)) +
  geom_errorbar(aes(ymin=0.01, ymax=Average+SD),
                width=.08, position=position_dodge()) +
  geom_bar(stat="identity") + scale_y_continuous(limits = c(0, 1.0)) + theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), legend.position = "none") 
plot(plot)
