library(Seurat)
library(ggplot2)
library(here)

# Input
cov02 <- Read10X("~/work/Data/cov18")

# Analysis
all.equal(colnames(cov02[["Gene Expression"]]), colnames(cov02[["Antibody Capture"]]))

cov02_comb <- CreateSeuratObject(counts = cov02[["Gene Expression"]])
cov02_comb[["ADT"]] <- CreateAssayObject(counts = cov02[["Antibody Capture"]])

Assays(cov02_comb)

# First working with RNA data
DefaultAssay(cov02_comb) <- "RNA"
cov02_comb[["percent.mt"]] <- PercentageFeatureSet(cov02_comb, pattern = "^MT-")
violin_RNA_QC <- VlnPlot(cov02_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, assay = "RNA")

cov02_comb <- subset(cov02_comb, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

cov02_comb <- NormalizeData(cov02_comb)
cov02_comb <- FindVariableFeatures(cov02_comb)
cov02_comb <- ScaleData(cov02_comb)
cov02_comb <- RunPCA(cov02_comb, verbose = FALSE)
cov02_comb <- FindNeighbors(cov02_comb, dims = 1:30)
cov02_comb <- FindClusters(cov02_comb, resolution = 1, verbose = FALSE)
cov02_comb <- RunUMAP(cov02_comb, dims = 1:30)
UMAP_RNA0.5 <- DimPlot(cov02_comb, label = TRUE, pt.size = 0.01)

# Output
saveRDS(UMAP_RNA0.5, file = "~/work/figures/UMAP_RNA0.5.rds")
saveRDS(violin_RNA_QC, file = "~/work/figures/violin_RNA_QC.rds")
ggsave(file.path("~/work/figures","UMAP_RNA0.5.png"),plot = UMAP_RNA0.5)
saveRDS(cov02_comb , file = "~/work/Data/Cov2RNA_clustered.rds")

