library(vegan)
library(tidyverse)
library(ggplot2)
library(DESeq2)

#otu_table input
otu <- read.csv("ASV_table.csv",row.names = 1, header = T)

#set threshold to remove rare species
#briefly check rowSums
mean(rowSums(otu))
max(rowSums(otu))
min(rowSums(otu))
otu = otu[which(rowSums(otu) > 100),]

#the row name of otu_table need to be sample id
otu <- t(otu)

#group info input
#the group file can be output using the kmeans scripts
group_info <- read.csv("cluster_member.csv", header=T, row.names = 1)
group_info <- group_info %>%
  mutate(sample_id = row.names(group_info))

#make otu_sub for analysis
#the definition of control and experimental groups does not impact SIMPER,
# but will impact the following differential abundance calculation
group_sub <- group_info %>%
  filter(cluster_id == 1 | cluster_id == 2) %>%
  filter(str_detect(sample_id,"RZ") == FALSE) %>% #<------ choose RT or RZ (Invert selection)
  mutate(new_cluster_id = case_when(cluster_id == 2  ~ "con", # define control group
                                    cluster_id == 1  ~ "exp")) # define experimental group
group_sub$new_cluster_id <- as.factor(group_sub$new_cluster_id)

otu_sub <- otu[which(rownames(otu)%in%group_sub$sample_id),]

#rarefy otu_sub
#this step is a MUST DO step!
#check reads number of each sample
rowSums(otu_sub)
mean(rowSums(otu_sub))
min(rowSums(otu_sub))
#remove samples with low reads if needed
otu_sub = otu_sub[which(rowSums(otu_sub) > 10000),]
rowSums(otu_sub)
#keep group_sub and otu_sub have same samples
group_sub <- group_sub[which(group_sub$sample_id%in%row.names(otu_sub)),]

#set rarefy number
min_read_num <- min(rowSums(otu_sub))
#rarefy otu
otu_sub = rrarefy(otu_sub, min_read_num)
rowSums(otu_sub)
#(optional) check the rarecurve
#rarecurve(otu_sub, col = "grey")

#Applying SIMPER method
#https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/simper/
#Similarity percentage (SIMPER) partitions the Bray-Curtis dissimilarity for every pair of sample units,
# and then calculates the average contribution of each species to the difference between the sample units.
#average: Average contribution to overall dissimilarity
#Permutation p-value: Probability of getting a larger or equal average contribution 
# in random permutation of the group factor.
sim_res <- simper(otu_sub, group_sub$new_cluster_id)

#arrange results
sim_res_sum <- summary(sim_res)
sim_res_sum <- sim_res_sum$con_exp %>%
  select(average, p)
sim_res_sum <- sim_res_sum %>%
  mutate(otu_id = row.names(sim_res_sum))
sim_res_sum$average <- round(sim_res_sum$average,3)

#Applying DESeq2 to detect differential abundance species
#make a read counts table
deseq2_otu_data <- as.data.frame(t(otu_sub + 1))
#check if the columns of the count data in the same order as rows names of the sample mapping
all(colnames(deseq2_otu_data) == rownames(group_sub))
#make deseq2 matrix
dds <- DESeqDataSetFromMatrix(countData=deseq2_otu_data, colData=group_sub, design= ~ new_cluster_id)
#run deseq function
dds <- DESeq(dds)
#extract results
res <- data.frame(results(dds))
Diff_res <- res %>%
  select(baseMean,log2FoldChange,pvalue,padj) %>%
  mutate(otu_id = row.names(res))

#relative abundance data
#only keep samples from experimental group
asv_table_reads <- as.data.frame(t(otu_sub))
asv_table_reads <- asv_table_reads[,group_sub$new_cluster_id == "exp"]

#normalize to percentage data
asv_table_colsum <- as.data.frame(colSums(asv_table_reads)) 
asv_table_colsum <- as.data.frame(t(asv_table_colsum))

asv_table_percentage <- asv_table_reads
for (i in (1:length(asv_table_reads[,1]))) {
  asv_table_percentage[i,] <- asv_table_reads[i,]/asv_table_colsum[1,]
}
relative_abundance <- as.data.frame(rowMeans(asv_table_percentage))
names(relative_abundance) <- "relative_abund"
relative_abundance <- relative_abundance %>%
  mutate(otu_id = row.names(relative_abundance))

#merge and get full results
All_res <- inner_join(sim_res_sum, Diff_res, by='otu_id')
All_res <- inner_join(All_res, relative_abundance, by='otu_id')

All_res_plot_simp <- All_res %>%
  filter(log2FoldChange >=  0) %>%
  filter(relative_abund >  0) %>%
  arrange(desc(average))

#simple plot and check
simplot <- ggplot(data = All_res_plot_simp, 
                     mapping = aes(x = average, y = log2FoldChange, colour = relative_abund )) + 
  scale_color_gradient2(mid = "#99CCFF" , high = "red") +
  geom_point(aes(size = relative_abund), alpha = 5/10)

simplot

#custumize plot
All_res_plot <- All_res_plot_simp %>%
  mutate(log2FC_new = case_when(log2FoldChange >=  10  ~ 10,
                                log2FoldChange < 10  ~ log2FoldChange)) %>%
  mutate(rel_abu_new = case_when(relative_abund >=  0.1  ~ 0.1,
                                 relative_abund < 0.1  ~ relative_abund)) %>%
  mutate(padj_new = case_when(padj <=  0.001  ~ 0.001,
                              padj <=  0.01  ~ 0.01,
                              padj <=  0.1  ~ 0.1,
                              padj <=  1  ~ 1))

#manual ruler
ruler <- All_res_plot[1,]
ruler[1,1] <- 0.08 #average x
ruler[1,9] <- 10 #log2FC_new y
ruler[1,11] <- 0.001 #padj color
ruler[1,10] <- 0.0001 #rel_abu_new size

All_res_plot <- rbind(All_res_plot, ruler)

plot_res_2 <- ggplot(data = All_res_plot, 
                     mapping = aes(x = average, y = log2FC_new, colour = rel_abu_new )) + 
  scale_color_gradient2(mid = "#99CCFF" , high = "#CC3300") +
  geom_point(aes(size = rel_abu_new), alpha = 9/10) +
  scale_size(range = c(0.5, 5)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "#666666", size = 0.8) +
  geom_vline(xintercept = 0.02, linetype = "dashed", color = "#666666", size = 0.8) +
  #scale_x_continuous(limits = c(0,0.1)) +
  #scale_y_continuous(limits = c(0,13)) +
  theme_light()+
  theme(aspect.ratio = 1/1)
plot_res_2

write.csv(All_res_plot,"OU_All_res_plot.csv", quote = FALSE)
ggsave("OU_diff_abun_ASV.pdf", plot = plot_res_2, width = 5, height = 5, dpi=600)

