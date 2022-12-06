library(Seurat)
library(ggplot2)
library(dplyr)

# Input
cov02_comb <- readRDS(file = "~/work/Data/Cov2RNA_clustered.rds")

# Analysis

DefaultAssay(cov02_comb) <- "ADT"

ADT_hist_CD3 <- hist(cov02_comb@assays[["ADT"]]@counts["CD3--UCHT1-TSA", ], main = "Antibody count of gene CD3", xlim = c(0,800), breaks = 2000,
     xlab = "Count")
RNA_hist_CD3 <- hist(cov02_comb@assays[["RNA"]]@counts["CD3E", ], main = "RNA count of gene CD3E",
     xlab = "Count")


ADR_10sample_boxplot <- boxplot( as.matrix(cov02_comb@assays[["ADT"]]@counts[, sample(c(1:1399), 10,  replace = F)]),
        main = "Antibody count of 10 randomly drawn cells", xlab = "Cell", ylab = "Antibody count", xaxt="n")

# Output
saveRDS(ADT_hist_CD3 , file = "~/work/figures/ADT_hist_CD3.rds")
saveRDS(RNA_hist_CD3 , file = "~/work/figures/RNA_hist_CD3.rds")
saveRDS(ADR_10sample_boxplot, file = "~/work/figures/ADR_10sample_boxplot.rds")

