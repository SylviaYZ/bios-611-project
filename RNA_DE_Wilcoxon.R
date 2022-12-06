library(Seurat)
library(ggplot2)
library(dplyr)

# Input
cov02_comb <- readRDS(file = "~/work/Data/Cov2RNA_clustered.rds")

# Analysis
DE_cluster <- FindAllMarkers(cov02_comb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DE_top10 <- DE_cluster %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) 
RNA_DEcluster_heatmap <- DoHeatmap(cov02_comb, features = DE_top10$gene) + NoLegend() + ggtitle("Top 10 DE genes in each cluster of cells")

# Output
saveRDS(RNA_DEcluster_heatmap, file = "~/work/figures/RNA_DEcluster_heatmap.rds")
