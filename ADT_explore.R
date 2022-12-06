library(Seurat)
library(ggplot2)
library(dplyr)

# Input
cov02_comb <- readRDS(file = "~/work/Data/Cov2RNA_clustered.rds")

# Analysis

DefaultAssay(cov02_comb) <- "ADT"

ADT_hist_CD3 <- hist(cov02_comb@assays[["ADT"]]@counts["CD3--UCHT1-TSA", ], main = "Antibody count of gene CD3", xlim = c(0,800), breaks = 200,
     xlab = "Count")
RNA_hist_CD3 <- hist(cov02_comb@assays[["RNA"]]@counts["CD3E", ], main = "RNA count of gene CD3E",
     xlab = "Count")

# Output
saveRDS(ADT_hist_CD3 , file = "~/work/figures/ADT_hist_CD3.rds")
saveRDS(RNA_hist_CD3 , file = "~/work/figures/RNA_hist_CD3.rds")
saveRDS(as.matrix(cov02_comb@assays[["ADT"]]@counts[, sample(c(1:1399), 10,  replace = F)]),"~/work/Data/ADT_10sample.rds" )

