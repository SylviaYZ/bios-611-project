library(Seurat)
library(ggplot2)
library(dplyr)

# Input
cov02_comb <- readRDS(file = "~/work/Data/Cov2RNA_clustered.rds")

# Add label
cov02_comb$Celltype_RNA <- as.factor(case_when(cov02_comb$seurat_clusters %in% c(0, 3, 8) ~ "T cell",
                                 cov02_comb$seurat_clusters %in% c(1,2, 11) ~ "NK cell",
                                 cov02_comb$seurat_clusters %in% c(4,5,6) ~ "CD14 monocyte & DC",
                                 cov02_comb$seurat_clusters %in% c(7, 13) ~ "CD8 T cell",
                                 cov02_comb$seurat_clusters %in% c(9) ~ "CD4 T cell",
                                 cov02_comb$seurat_clusters %in% c(10) ~ "B cell",
                                 cov02_comb$seurat_clusters %in% c(11) ~ "Platelet",
                                 TRUE ~ "Unknown"))

UMAP_celltypeRNA <- UMAPPlot(cov02_comb, group.by = "Celltype_RNA", label = F, pt.size = 0.01) 

cov02_comb$Celltype_ADT <- as.factor(case_when(cov02_comb$seurat_clusters %in% c(3, 8) ~ "T cell",
                                               cov02_comb$seurat_clusters %in% c(1,2, 11) ~ "NK cell",
                                               cov02_comb$seurat_clusters %in% c(4) ~ "CD14 monocyte & DC",
                                               cov02_comb$seurat_clusters %in% c(5,6) ~ "CD14 monocyte & DC",
                                               cov02_comb$seurat_clusters %in% c(7, 13) ~ "CD8 T cell",
                                               cov02_comb$seurat_clusters %in% c(0, 8, 9) ~ "CD4 T cell",
                                               cov02_comb$seurat_clusters %in% c(10) ~ "B cell",
                                               cov02_comb$seurat_clusters %in% c(11) ~ "Platelet",
                                               TRUE ~ "Unknown"))

UMAP_celltypeADT <- UMAPPlot(cov02_comb, group.by = "Celltype_ADT", label = F, pt.size = 0.01) 



# Output
ggsave("~/work/UMAP_celltypeRNA.png" , plot = UMAP_celltypeRNA)
ggsave("~/work/UMAP_celltypeADT.png" , plot = UMAP_celltypeADT)
