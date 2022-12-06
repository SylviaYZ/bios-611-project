library(Seurat)
library(ggplot2)

# Input
cov02_comb <- readRDS(file = "~/work/Data/Cov2RNA_clustered.rds")

# Plots
RNA_violin_T <- VlnPlot(cov02_comb, features = c("CD3D","CD3E"))
RNA_violin_T4_T8 <- VlnPlot(cov02_comb, features = c("CD4","CD8A"))
RNA_violin_B <- VlnPlot(cov02_comb, features = c("MS4A1","CD19"))
RNA_violin_NK <- VlnPlot(cov02_comb, features = c("GNLY","NKG7"))
RNA_violin_CD14mono <- VlnPlot(cov02_comb, features = c("CD14","LYZ"))
RNA_violin_DC <- VlnPlot(cov02_comb, features = c("FCER1A","CST3"))
RNA_violin_PL <- VlnPlot(cov02_comb, features = c("PPBP"))

# Output
saveRDS(RNA_violin_T, file = "~/work/figures/RNA_violin_T.rds")
saveRDS(RNA_violin_T4_T8, file = "~/work/figures/RNA_violin_T4_T8.rds")
saveRDS(RNA_violin_B, file = "~/work/figures/RNA_violin_B.rds")
saveRDS(RNA_violin_NK, file = "~/work/figures/RNA_violin_NK.rds")
saveRDS(RNA_violin_CD14mono, file = "~/work/figures/RNA_violin_CD14mono.rds")
saveRDS(RNA_violin_DC, file = "~/work/figures/RNA_violin_DC.rds")
saveRDS(RNA_violin_PL, file = "~/work/figures/RNA_violin_PL.rds")