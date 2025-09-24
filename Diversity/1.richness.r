library(picante)
library(vegan)
library(tidyverse)

#define function for calculating alpha diversity
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}

#input asv_table from qiime2
otu <- read.csv("ASV_table.csv",row.names = 1, header = T)

otu <- t(otu)
#only keep rows with rowsum > 0
otu=subset(otu,rowSums(otu) > 0)

#rarefy otu
#check reads number of each sample
rowSums(otu)
mean(rowSums(otu))
min(rowSums(otu))
#remove samples with low reads if needed
otu = otu[which(rowSums(otu) > 20000),]
rowSums(otu)

#set rarefy number
min_read_num <- min(rowSums(otu))
#rarefy otu
otu = rrarefy(otu, min_read_num)
rowSums(otu)
#(optional) check the rarecurve
#rarecurve(otu, col = "#0000004D")

#calculating alpha diversity:
#Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage
alpha_all <- alpha(otu, base = 2)
alpha_all <- alpha_all %>%
  mutate(id = row.names(alpha_all))

#read metadata
metadata <- read.delim('sample-metadata.tsv', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
metadata <- metadata %>%
  mutate(id = row.names(metadata))

alpha_all <- inner_join(alpha_all, metadata)

#output results
write.csv(alpha_all, 'OU_alpha.csv', quote = FALSE)
