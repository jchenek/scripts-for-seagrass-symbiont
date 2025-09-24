library(ggplot2)
library(reshape2)
library(tidyverse)

#read asv table
asv_table_reads <- read.csv("ASV_table.csv",row.names = 1, header = T)
asv_table_reads <- asv_table_reads %>%
  mutate(OTU = row.names(asv_table_reads))

#read target asvs
target_asvs <- read.delim('target_asvs.txt', row.names = NULL, sep = '\t', 
                          stringsAsFactors = FALSE, check.names = FALSE, header = F)
names(target_asvs) <- c("OTU","taxo")

#create target asv table
target_asv_table_reads <- inner_join(asv_table_reads,target_asvs)

#make pure asv table
row.names(asv_table_reads) <- asv_table_reads$OTU
asv_table_reads <- asv_table_reads %>%
  select(-OTU)
row.names(target_asv_table_reads) <- target_asv_table_reads$taxo
target_asv_table_reads <- target_asv_table_reads %>%
  select(-OTU,-taxo)

#remove Singleton 
asv_table_reads <- asv_table_reads[which(rowMeans(asv_table_reads)> 1) ,]

#get total reads of each sample
asv_table_colsum <- as.data.frame(colSums(asv_table_reads)) 
asv_table_colsum <- as.data.frame(t(asv_table_colsum))

#normalize all to percentage data
asv_table_percentage <- asv_table_reads
for (i in (1:length(asv_table_reads[,1]))) {
  asv_table_percentage[i,] <- asv_table_reads[i,]/asv_table_colsum[1,]
}

#normalize target to percentage data
target_asv_table_reads_percentage <- target_asv_table_reads
for (i in (1:length(target_asv_table_reads[,1]))) {
  target_asv_table_reads_percentage[i,] <- target_asv_table_reads[i,]/asv_table_colsum[1,]
}

#check target percentage
colSums(target_asv_table_reads_percentage)

#check top 10 asvs
#RT
asv_table_percentage_RT <- asv_table_percentage %>%
  select(matches("RT"))
asv_table_percentage_RT_mean <- as.data.frame(rowMeans(asv_table_percentage_RT))
asv_table_percentage_RT_mean <- asv_table_percentage_RT_mean %>%
 arrange(desc(`rowMeans(asv_table_percentage_RT)`)) 
#write.csv(asv_table_percentage_RT_mean, "OU_asv_table_percentage_RT_mean.csv")
#RZ
asv_table_percentage_RZ <- asv_table_percentage %>%
  select(matches("RZ"))
asv_table_percentage_RZ_mean <- as.data.frame(rowMeans(asv_table_percentage_RZ))
asv_table_percentage_RZ_mean <- asv_table_percentage_RZ_mean %>%
  arrange(desc(`rowMeans(asv_table_percentage_RZ)`))
#write.csv(asv_table_percentage_RZ_mean, "OU_asv_table_percentage_RZ_mean.csv")
#write.csv(target_asv_table_reads_percentage, "OU_target_asv_table_reads_percentage.csv")

#make target plot
#for RT
target_asv_table_reads_percentage_RT <- target_asv_table_reads_percentage %>%
  select(matches("RT"))
last_row_RT <- as.data.frame(t(1-colSums(target_asv_table_reads_percentage_RT)))
row.names(last_row_RT) <- "AAothers"
target_asv_table_reads_percentage_RT <- rbind(last_row_RT,target_asv_table_reads_percentage_RT)
target_asv_table_reads_percentage_RT <- target_asv_table_reads_percentage_RT %>%
  mutate(level_id = row.names(target_asv_table_reads_percentage_RT))
  
gg_data_r <- melt(target_asv_table_reads_percentage_RT, id.vars="level_id", variable.name="sample_id", value.name="abundance")
#simple plot and check
# ggplot(gg_data_r, aes(sample_id, abundance, fill=level_id)) +
#   geom_bar(stat="identity",colour = 'grey',linewidth = 0.05) +
#   guides(fill=guide_legend(reverse=F)) +
#   theme(panel.grid = element_blank()) +
#   theme_bw()
#design color and plot
#get color code
#adjustcolor("grey", alpha.f = 0.1)
gg_plot_r <- ggplot(gg_data_r, aes(sample_id, abundance, fill=level_id)) +
  geom_bar(stat="identity",colour = '#BEBEBE4D',linewidth = 0.1) +
  guides(fill=guide_legend(reverse=F)) +
  scale_fill_manual(values =  rev(c(
    '#9999CC','#AECDE0','#3C76AF',
    '#F0B770','#70BAA8','#E08985','#BEBEBE1A'))) +
  theme(panel.grid = element_blank()) + 
  theme_void()
gg_plot_r
ggsave("OU_taxo_barplot_RT.pdf", plot = gg_plot_r, width = 2.5, height = 4, dpi=600)

#check mean relative abundance
gb_r <- gg_data_r %>%
  group_by(level_id) %>%
  summarise(meanS = mean(abundance))

#for RZ
target_asv_table_reads_percentage_RZ <- target_asv_table_reads_percentage %>%
  select(matches("RZ"))
last_row_RZ <- as.data.frame(t(1-colSums(target_asv_table_reads_percentage_RZ)))
row.names(last_row_RZ) <- "AAothers"
target_asv_table_reads_percentage_RZ <- rbind(last_row_RZ,target_asv_table_reads_percentage_RZ)
target_asv_table_reads_percentage_RZ <- target_asv_table_reads_percentage_RZ %>%
  mutate(level_id = row.names(target_asv_table_reads_percentage_RZ))

gg_data_s <- melt(target_asv_table_reads_percentage_RZ, id.vars="level_id", variable.name="sample_id", value.name="abundance")
#simple plot and check
# ggplot(gg_data_s, aes(sample_id, abundance, fill=level_id)) +
#   geom_bar(stat="identity",colour = 'grey',linewidth = 0.05) +
#   guides(fill=guide_legend(reverse=F)) +
#   theme(panel.grid = element_blank()) +
#   theme_bw()
#design color and plot
#get color code
#adjustcolor("grey", alpha.f = 0.1)
gg_plot_s <- ggplot(gg_data_s, aes(sample_id, abundance, fill=level_id)) +
  geom_bar(stat="identity",colour = '#BEBEBE4D',linewidth = 0.1) +
  guides(fill=guide_legend(reverse=F)) +
  scale_fill_manual(values =  rev(c(
    '#9999CC','#AECDE0','#3C76AF',
    '#F0B770','#70BAA8','#E08985','#BEBEBE1A'))) +
  theme(panel.grid = element_blank()) + 
  theme_void()
gg_plot_s
ggsave("OU_taxo_barplot_RZ.pdf", plot = gg_plot_s, width = 2.5, height = 4, dpi=600)

#check mean relative abundance
gb_s <- gg_data_s %>%
  group_by(level_id) %>%
  summarise(meanS = mean(abundance))

#for control
target_asv_table_reads_percentage_cont <- target_asv_table_reads_percentage %>%
  select(matches("NS"),matches("BS"),matches("RS")) 

last_row_cont <- as.data.frame(t(1-colSums(target_asv_table_reads_percentage_cont)))
row.names(last_row_cont) <- "AAothers"
target_asv_table_reads_percentage_cont <- rbind(last_row_cont,target_asv_table_reads_percentage_cont)
target_asv_table_reads_percentage_cont <- target_asv_table_reads_percentage_cont %>%
  mutate(level_id = row.names(target_asv_table_reads_percentage_cont))

gg_data_c <- melt(target_asv_table_reads_percentage_cont, id.vars="level_id", variable.name="sample_id", value.name="abundance")
#simple plot and check
# ggplot(gg_data_c, aes(sample_id, abundance, fill=level_id)) +
#   geom_bar(stat="identity",colour = 'grey',linewidth = 0.05) +
#   guides(fill=guide_legend(reverse=F)) +
#   theme(panel.grid = element_blank()) +
#   theme_bw()
#design color and plot
#get color code
#adjustcolor("grey", alpha.f = 0.1)
gg_plot_c <- ggplot(gg_data_c, aes(sample_id, abundance, fill=level_id)) +
  geom_bar(stat="identity",colour = '#BEBEBE4D',linewidth = 0.1) +
  guides(fill=guide_legend(reverse=F)) +
  scale_fill_manual(values =  rev(c(
    '#9999CC','#AECDE0','#3C76AF',
    '#F0B770','#70BAA8','#E08985','#BEBEBE1A'))) +
  theme(panel.grid = element_blank()) + 
  theme_void()
gg_plot_c
ggsave("OU_taxo_barplot_Sed.pdf", plot = gg_plot_c, width = 2.5, height = 4, dpi=600)

#check mean relative abundance
gb_c <- gg_data_c %>%
  group_by(level_id) %>%
  summarise(meanS = mean(abundance))
