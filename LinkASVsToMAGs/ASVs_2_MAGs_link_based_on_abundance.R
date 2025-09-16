library(ggplot2)
library(reshape2)
library(tidyverse)
library(psych)
library(corrplot)

#This scripts is trying to figure out the source MAG of each ASV
# based on the correlation between ASVs-abund and MAGs-abund in many samples 
#In the end of this scripts
#check OU_spearman_res_r.csv, the most correlated MAG-ASV should be coupled
#check OU_lm_plots.pdf, the most correlated (should also be the only highlighted) MAG-ASV should be coupled
#these results should be consist with each other

#read MAGs abund data from coverm 
coverm_res <- read.delim('IN_coverm.tsv', row.names = NULL, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
row.names(coverm_res) <- coverm_res$Genome
coverm_res <- coverm_res %>%
  select(-Genome)

coverm_res <- as.data.frame(t(coverm_res))
coverm_res <- coverm_res %>%
  select(-unmapped) %>%
  mutate(bam_file = row.names(coverm_res))

#input total mapped reads number
coverm_reads <- read.delim('IN_bam_total_reads', row.names = NULL, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#only keep samples > 0.1M
coverm_reads <- coverm_reads %>%
  filter(mapped_reads > 100000)

#join coverm data
coverm_res_num <- inner_join(coverm_reads, coverm_res)
row.names(coverm_res_num) <- coverm_res_num$bam_file
coverm_res_num <- coverm_res_num %>%
  select(-bam_file, -mapped_reads)
coverm_res_num <- as.data.frame(t(coverm_res_num))

#normalize to percentage data
input_data_colsum <- as.data.frame(colSums(coverm_res_num)) 
input_data_colsum <- as.data.frame(t(input_data_colsum))

ratio_data <- coverm_res_num
for (i in (1:length(coverm_res_num[,1]))) {
  ratio_data[i,] <- coverm_res_num[i,]/input_data_colsum[1,]
}
write.csv(ratio_data, "OU_ratio_data_full.csv")

#read metadata
metadata <- read.delim('IN_sample-metadata.tsv', row.names = NULL, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#only keep RT and RZ samples
metadata_RT_RZ <- metadata %>%
  filter(type == "RT" | type == "RZ") %>%
  arrange(desc(type))

metadata_RT_RZ_id <- as.data.frame(metadata_RT_RZ$`sample-id`)
names(metadata_RT_RZ_id) = "sample_id"

#get RT RZ ratio
ratio_data_RT_RZ <- as.data.frame(t(ratio_data))
ratio_data_RT_RZ <- ratio_data_RT_RZ %>%
  mutate(sample_id = row.names(ratio_data_RT_RZ))

ratio_data_RT_RZ <- inner_join(metadata_RT_RZ_id,ratio_data_RT_RZ )
row.names(ratio_data_RT_RZ) <- ratio_data_RT_RZ$sample_id
ratio_data_RT_RZ <- ratio_data_RT_RZ %>%
  select(-sample_id)
ratio_data_RT_RZ <- as.data.frame(t(ratio_data_RT_RZ))

write.csv(ratio_data_RT_RZ, "OU_ratio_data_RT_RZ.csv")

#read target MAGs
MAGs_id <- read.delim('IN_MAGs_id.tsv', row.names = NULL, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#only keep target MAGs
ratio_data_RT_RZ_MAGs <- ratio_data_RT_RZ %>%
  mutate(MAG_id = row.names(ratio_data_RT_RZ))

ratio_data_RT_RZ_MAGs <- inner_join(MAGs_id, ratio_data_RT_RZ_MAGs)
row.names(ratio_data_RT_RZ_MAGs) <- ratio_data_RT_RZ_MAGs$MAG_id
ratio_data_RT_RZ_MAGs <- ratio_data_RT_RZ_MAGs %>%
  select(-MAG_id)

write.csv(ratio_data_RT_RZ_MAGs, "OU_ratio_data_RT_RZ_MAGs.csv")

#input asv table
asv_table <- read.csv('IN_clean_asv_table.csv', row.names = NULL)
row.names(asv_table) <- asv_table$OTU
asv_table <- asv_table %>%
  select(-OTU)
#record sample total reads
asv_table_colsum <- as.data.frame(colSums(asv_table)) 
asv_table_colsum <- as.data.frame(t(asv_table_colsum))

#input target asv
ASVs_id <- read.delim('IN_ASVs_id.tsv', row.names = NULL, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#extract target asv table
asv_table <- read.csv('IN_clean_asv_table.csv', row.names = NULL)
ASVs_table <- inner_join(ASVs_id, asv_table)
row.names(ASVs_table) <- ASVs_table$OTU
ASVs_table_reads <- ASVs_table %>%
  select(-OTU)

#normalize to percentage
ASVs_table <- ASVs_table_reads
for (i in (1:length(ASVs_table_reads[,1]))) {
  ASVs_table[i,] <- ASVs_table_reads[i,]/asv_table_colsum[1,]
}

#inner join asv and mag data
ratio_data_RT_RZ_MAGs <- as.data.frame(t(ratio_data_RT_RZ_MAGs))
ratio_data_RT_RZ_MAGs_samples <- as.data.frame(row.names(ratio_data_RT_RZ_MAGs)) 
names(ratio_data_RT_RZ_MAGs_samples) <- "samples"

ASVs_table <- as.data.frame(t(ASVs_table))
ASVs_table <- ASVs_table %>%
  mutate(samples = row.names(ASVs_table))

ASVs_table = inner_join(ratio_data_RT_RZ_MAGs_samples, ASVs_table)
row.names(ASVs_table) <- ASVs_table$samples
ASVs_table = ASVs_table %>%
  select(-samples)

#check the valid samples used for calculating correlation, this number should be big enough (>15)
length(ASVs_table[,1])

#get spearman table
spearman_res_r <- data.frame()
spearman_res_p <- data.frame()

for (i in (1:(length(ASVs_table[1,])))) {
  for(j in (1:length(ratio_data_RT_RZ_MAGs[1,]))){
    corr_matrix <- corr.test(ASVs_table[,i], ratio_data_RT_RZ_MAGs[,j], method = 'spearman')
    spearman_res_r[j,i] <- round(corr_matrix$r,2)
    spearman_res_p[j,i] <- round(corr_matrix$p,2)
  }
}

names(spearman_res_r) <- names(ASVs_table)
names(spearman_res_p) <- names(ASVs_table)
row.names(spearman_res_r) <- names(ratio_data_RT_RZ_MAGs)
row.names(spearman_res_p) <- names(ratio_data_RT_RZ_MAGs)

#plot cor
spearman_res_r_plot <- as.matrix(t(-spearman_res_r))
corrplot(spearman_res_r_plot,method = 'pie',
                          addgrid.col = "#BEBEBE4D",
                          tl.pos = 'n',
                          cl.pos = 'n'
                          )
#M <- cor(mtcars)
#corrplot(M)

#check OU_spearman_res_r.csv, the most correlated MAG-ASV should be coupled
write.csv(spearman_res_r, "OU_spearman_res_r.csv")
write.csv(spearman_res_p, "OU_spearman_res_p.csv")

#check OU_lm_plots.pdf, the most correlated (should also be the only highlighted)
# MAG-ASV should be coupled
#regression plot
ncol_plot = length(ratio_data_RT_RZ_MAGs[1,])
nrow_plot = length(ASVs_table[1,])

pdf(file = "OU_lm_plots.pdf", wi = 60, he = 11)
par(mfrow = c(nrow_plot,ncol_plot),mar=c(1,0.5,1,0.5),oma=c(1,1,1,1))
i = 1
while (i <= nrow_plot) {
  j = 1
  while(j <= ncol_plot) {
    ASV_name <- names(ASVs_table)[i]
    MAG_name <- names(ratio_data_RT_RZ_MAGs)[j]
    fit_data <- as.data.frame(cbind(ASVs_table[,i],ratio_data_RT_RZ_MAGs[,j]))
    fit <- lm( V1 ~ V2, data = fit_data)
    #get R squared and p value
    print(summary(fit))
    #cancel scientific notation
    options(scipen=200)
    r_sqa <- round(summary(fit)$r.squared,2)
    p_val <- round(summary(fit)$coefficients[2,4],4)
    #plot
    plot(x=fit_data$V2, y=fit_data$V1,
         col = "#6699FF80", pch= 16, cex=3.5,
         #main = paste(MAG_name,"\n",ASV_name,sep = ""),cex.main=0.5,
         xlab = "", ylab= "",
         xaxt = "n", yaxt = "n")
    if(p_val <= 0.01 & r_sqa >= 0.7){
      abline(fit, col = "#FF33334D", lty = 1, lwd = 20, )
    }
    j = j + 1
  }
  i = i + 1  
}
dev.off()
