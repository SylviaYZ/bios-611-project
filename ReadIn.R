library(Seurat)
library(ggplot2)
library(here)
cov01 <- Read10X("~/work/Data/cov01")

all.equal(colnames(cov01[["Gene Expression"]]), colnames(cov01[["Antibody Capture"]]))

cov01_comb <- CreateSeuratObject(counts = cov01[["Gene Expression"]])
cov01_comb[["ADT"]] <- CreateAssayObject(counts = cov01[["Antibody Capture"]])

Assays(cov01_comb)

# First working with RNA data
DefaultAssay(cov01_comb) <- "RNA"
cov01_comb <- NormalizeData(cov01_comb)
cov01_comb <- FindVariableFeatures(cov01_comb)
cov01_comb <- ScaleData(cov01_comb)
cov01_comb <- RunPCA(cov01_comb, verbose = FALSE)
cov01_comb <- FindNeighbors(cov01_comb, dims = 1:30)
cov01_comb <- FindClusters(cov01_comb, resolution = 0.8, verbose = FALSE)
cov01_comb <- RunUMAP(cov01_comb, dims = 1:30)
UMAP_RNA0.8 <- DimPlot(cov01_comb, label = TRUE)
saveRDS(UMAP_RNA0.8, file = "~/work/figures/UMAP_RNA0.8.rds")
ggsave(file.path("~/work/figures","UMAP_RNA0.8.png"),plot = UMAP_RNA0.8)

# Now working with RNA data
cov01_comb <- NormalizeData(cov01_comb, normalization.method = "CLR", margin = 2, assay = "ADT")

DefaultAssay(cov01_comb) <- "ADT"
p1 <- FeaturePlot(cov01_comb, "CD3--UCHT1-TSA", cols = c("lightgrey", "darkgreen")) + ggtitle("CD3 protein")
DefaultAssay(cov01_comb) <- "RNA"
p2 <- FeaturePlot(cov01_comb, "CD3") + ggtitle("CD3 RNA")

p1 | p2

adt_markers <- FindAllMarkers(cov01_comb, assay = "ADT")
rna_markers <- FindAllMarkers(cov01_comb, assay = "RNA")
