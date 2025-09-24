####################################################################################################################
library(tidyverse)
library(ggplot2)

#input diversity data
joined_data <- read.csv("OU_alpha.csv", header=T)
joined_data <- joined_data %>%
  select(Richness,type)

#plot
#simple plot
ggplot(joined_data, aes(x=type, y=Richness)) + 
  geom_boxplot() +
  theme(axis.text.x=element_text(angle=45,size=8))

#plot with dots
ggplot(joined_data, aes(x=type, y=Richness)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter()) +
  theme(axis.text.x=element_text(angle=45,size=8))

#plot with order
joined_data$type <- factor(joined_data$type, levels = c("NS","BS","RS","RT","RZ"))

p <- ggplot(joined_data, aes(x=type, y=Richness))  + 
  stat_boxplot(geom="errorbar",width=0.38,color="#333333",linewidth=1)+
  geom_boxplot(size = 1,color="#333333") +
  geom_jitter(alpha = 0.7, colour = "#666666", position=position_jitter(0.15),size=3) +
  theme_bw()
p

ggsave("OU_boxplot.pdf", plot = p, width = 5, height = 5, dpi=600)

