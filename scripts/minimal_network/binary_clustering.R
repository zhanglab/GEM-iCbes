raw_data <- read.table('../../../minimal_network/ethanol_strain_new_test/E1M-mn_50_combine_all_repeats.tsv',
                   sep = "\t",
                   row.names = 1)
colnames(raw_data) <- 1:1000

require(pheatmap)
require(RColorBrewer)

# ## take a look of 1000 simulations by the heatmap
# pheatmap(
#   data,
#   clustering_distance_rows = "binary",
#   clustering_distance_cols = "binary",
#   clustering_method = "ward.D2"
# )

# remove the core-essential and non-essential, only care the flexible genes
index <- rowSums(raw_data) > 0 & rowSums(raw_data) < ncol(raw_data) 
data <- raw_data[index, ]

# ## take a look of the the binary clustered heatmap
# data.heatmap <- pheatmap(
#   data,
#   clustering_distance_rows = "binary",
#   clustering_distance_cols = "binary",
#   clustering_method = "ward.D2"
# )


# # Install
# install.packages("FactoMineR")
# 
# # Load
# library("FactoMineR")

# install.packages("factoextra")

require(factoextra)
# calculate distence among 1000 MNs
data.dist <- dist(t(data), method = "binary")
write.table(1-as.matrix(data.dist), 
            "../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/1000MNs_similarity.tsv",
            sep = "\t",
            col.names = NA)
pdf('../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/MNs_distence_heatmap.pdf', height=16, width=16)
gene.heatmap <- pheatmap(1-as.matrix(data.dist),
                         clustering_distance_rows = data.dist,
                         clustering_distance_cols = data.dist,
                         clustering_method = "ward.D2")
dev.off()

# calculate distence of all flexible genes
gene.dist <- dist(data, method = "binary")
write.table(1-as.matrix(gene.dist), 
            "../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/gene_similarity.tsv",
            sep = "\t",
            col.names = NA)
pdf('../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/gene_distence_heatmap.pdf', height=16, width=16)
gene.heatmap <- pheatmap(1-as.matrix(gene.dist),
                         clustering_distance_rows = gene.dist,
                         clustering_distance_cols = gene.dist,
                         clustering_method = "ward.D2")
dev.off()


# [Method 1]: silhouette, maximize the average silhouette, to show the difference among different groups
pdf('../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/silhouette.pdf')
fviz_nbclust(t(data), FUNcluster = hcut, diss = data.dist)
dev.off()

# [Method 2]: Elbow, changing point of the total within-cluster sum of square, to show the difference within each group
pdf('../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/Elbow.pdf')
fviz_nbclust(t(data),
             FUNcluster = hcut,
             diss = data.dist,
             method = "wss")
dev.off()

# [Method 3]: gap statistics, compare the group result to null model, maximize the differences
set.seed(123)
pdf('../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/gap_statistic.pdf')
fviz_nbclust(
  t(data),
  FUNcluster = hcut,
  diss = data.dist,
  method = "gap_stat",
  nboot = 100
)
dev.off()


# choose a k number and cut them into k groups
data.hcut <- hcut(data.dist, 10)
pdf('../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/Clustering_k7.pdf')
fviz_dend(data.hcut, rect = T)
dev.off()

# make a heatmap with clustering at k groups
pdf('../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/Clustering_heatmap_k7.pdf', height=16, width=6)
data.heatmap <- pheatmap(
  data,
  clustering_distance_rows = "binary",
  clustering_distance_cols = "binary",
  clustering_method = "ward.D2",
  cutree_cols = 10,
  show_rownames = T,
  show_colnames = F,
  legend = F
)
dev.off()


write.table(data.hcut$cluster,
            "../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/cluster.txt",
            sep = "\t",
            row.names = 1:1000,
            col.names = F)


# count gene number of each MN in the same cluster
gene.count <- colSums(raw_data)
min(gene.count)
max(gene.count)
names(gene.count) <- paste0("E1M-mn_50_NOtfba_",1:1000,"_out.tsv")

gene.count[gene.count == min(gene.count)]

# group1 <- gene.count[data.hcut$cluster==1]
# group2 <- gene.count[data.hcut$cluster==2]
# group3 <- gene.count[data.hcut$cluster==3]
# group4 <- gene.count[data.hcut$cluster==4]
# group5 <- gene.count[data.hcut$cluster==5]
# 
# group5[group5==min(gene.count)]


pdf('../../../minimal_network/ethanol_strain_new_test/cluster_E1M-mn_50/boxplot_cluster.pdf', height=6, width=8)
boxplot(gene.count~data.hcut$cluster,
        xlab="Groups",
        ylab="Gene number in each minimal network",
        ylim=c(250, 280))
dev.off()





