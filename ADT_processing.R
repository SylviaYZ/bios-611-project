library(Seurat)
library(ggplot2)
library(dplyr)

# Input
cov02_comb <- readRDS(file = "~/work/Data/Cov2RNA_clustered.rds")

# TCGGGCAGTAAGCGGT-1 
# 1171 

# Analysis
cov02_comb <- NormalizeData(cov02_comb, normalization.method = "CLR", margin = 2, assay = "ADT")

saveRDS(as.matrix(cov02_comb@assays[["ADT"]]@data[, sample(c(1:1399), 10,  replace = F)]),"~/work/Data/ADT_10sample_normalzied.rds" )

# T cell
VlnPlot(cov02_comb, c("CD3--UCHT1-TSA"),  assay = "ADT")
VlnPlot(cov02_comb, c("CD8--RPA-T8-TSA", "CD4--OKT4-TSA"))

# B cell
VlnPlot(cov02_comb, c("CD19--HIB19-TSA"))


# Output
DefaultAssay(cov02_comb ) <- "ADT"

