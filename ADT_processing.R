library(Seurat)
library(ggplot2)
library(dplyr)

# Input
cov02_comb <- readRDS(file = "~/work/Data/Cov2RNA_clustered.rds")

# TCGGGCAGTAAGCGGT-1 
# 1171 

# Analysis
cov02_comb <- NormalizeData(cov02_comb, normalization.method = "CLR", margin = 2, assay = "ADT")

# Output data!!
saveRDS(as.matrix(cov02_comb@assays[["ADT"]]@data[, sample(c(1:1399), 10,  replace = F)]),"~/work/Data/ADT_10sample_normalzied.rds" )

### Violin plots
ADT_violin_T <- VlnPlot(cov02_comb, c("CD3--UCHT1-TSA", "CD28--CD28-2-TSA"),  assay = "ADT")
ADT_violin_T4_T8 <- VlnPlot(cov02_comb, c("CD8--RPA-T8-TSA", "CD4--OKT4-TSA"), assay = "ADT")
ADT_violin_T4naive <- VlnPlot(cov02_comb, c("CD45RA--HI100-TSA"), assay = "ADT")
ADT_violin_T4memeff <- VlnPlot(cov02_comb, c("CD127--A019D5-TSA"), assay = "ADT")
ADT_violin_B <- VlnPlot(cov02_comb, c("CD19--HIB19-TSA", "CD20--2H7-TSA"),  assay = "ADT")
ADT_violin_NK <- VlnPlot(cov02_comb, features = c("CD56--5-1H11-TSA"),  assay = "ADT")
ADT_violin_CD14mono <- VlnPlot(cov02_comb, features = c("CD14--M5E2-TSA"),  assay = "ADT")
ADT_violin_DC <- VlnPlot(cov02_comb, features = c("FCER1a--AER-37-TSA", "CD11c--S-HCL-3-TSA"),  assay = "ADT")


# Output
saveRDS(ADT_violin_T, file = "~/work/figures/ADT_violin_T.rds")
saveRDS(ADT_violin_T4_T8, file = "~/work/figures/ADT_violin_T4_T8.rds")
saveRDS(ADT_violin_T4naive, file = "~/work/figures/ADT_violin_T4naive.rds")
saveRDS(ADT_violin_T4memeff, file = "~/work/figures/ADT_violin_T4memeff.rds")
saveRDS(ADT_violin_B, file = "~/work/figures/ADT_violin_B.rds")
saveRDS(ADT_violin_NK, file = "~/work/figures/ADT_violin_NK.rds")
saveRDS(ADT_violin_CD14mono, file = "~/work/figures/ADT_violin_CD14mono.rds")
saveRDS(ADT_violin_DC, file = "~/work/figures/ADT_violin_DC.rds")
