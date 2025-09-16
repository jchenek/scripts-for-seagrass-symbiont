library(vegan)
library(factoextra)
library(tidyverse)

#set random number and read data
getwd()
set.seed(520)

#col: species; row: samples
#hclust will group the rows
otu <- read.csv("ASV_table.csv",row.names = 1, header = T)

#make cluster for col or row
otu <- t(otu)

#normalization data, various method available.
otu <- decostand(otu,method = "hellinger")
#calculating distance, various method available.
DisTan <- vegdist(otu, method="bray") 

#obtain best cluster number
#there are many methods to get the number: https://www.cnblogs.com/think90/p/7133753.html
#https://medium.com/analytics-vidhya/how-to-determine-the-optimal-k-for-k-means-708505d204eb
#There are three recommende methods
#Set the max k number
MaxKn = 5

#The Elbow Method
#Calculate the Within-Cluster-Sum of Squared Errors (WSS) for different values of k, 
# and choose the k for which WSS becomes first starts to diminish.
# In the plot of WSS-versus-k, this is visible as an elbow.
p_wss <- fviz_nbclust(otu, kmeans, k.max = MaxKn, method  = "wss", diss = DisTan)
p_wss
ggsave("OU_p_wss.pdf", plot = p_wss, width = 3, height = 5, dpi=600)

#The Silhouette Method
#The silhouette value measures how similar a point is to its own cluster (cohesion) compared to other clusters (separation).
# The range of the Silhouette value is between +1 and -1. 
# A high value is desirable and indicates that the point is placed in the correct cluster.
p_silhouette <- fviz_nbclust(otu, kmeans, k.max = MaxKn, method  = "silhouette", diss = DisTan)
p_silhouette
ggsave("OU_p_silhouette.pdf", plot = p_silhouette, width = 3, height = 5, dpi=600)

#the Gap-Statistics Method
#this method is much slower than the other two
#https://towardsdatascience.com/k-means-clustering-and-the-gap-statistics-4c5d414acd29
p_gap <- fviz_nbclust(otu, kmeans, k.max = MaxKn, method  = "gap_stat", diss = DisTan)
p_gap
ggsave("OU_p_gap.pdf", plot = p_gap, width = 3, height = 5, dpi=600)

#check members of each group
#Optimal K number
OptK = 2
#the default algorithm, by Hartigan and Wong (1979) may be the smartest algorithm
#nstart: if centers is a number, how many random sets should be chosen?
kmeans_res <- kmeans(otu, centers = OptK, nstart = 100,
                     algorithm = "Hartigan-Wong")

#check cluster_res if the clustering is reasonable
cluster_res <-as.data.frame(kmeans_res$cluster)
colnames(cluster_res) <- c("cluster_id")
write.csv(cluster_res, "OU_cluster_member.csv")

#Dendrogram plot
dend_res <- hclust(DisTan, "ward.D2")
p_dend <- fviz_dend(dend_res, k = OptK, lwd = 0.8,
          k_colors = c("#C2B2CD", "#45AB84"))
p_dend
ggsave("OU_p_dend.pdf", plot = p_dend, width = 4, height = 5, dpi=600)

#cluster plot
fviz_cluster(kmeans_res, data = otu,
             geom = c("point", "text"),
             show.clust.cent = T)

#Analysis of similarities (ANOSIM) 
Group <- kmeans_res$cluster
data.anosim = anosim(DisTan, grouping = Group, permutations = 999)
summary(data.anosim)
