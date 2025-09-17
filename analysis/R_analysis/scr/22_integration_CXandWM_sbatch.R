# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)

# read in the data --------------------------------------------------------
# load the LUT
LUT <- read_csv("../../data/LUT_sample_WM_CX.csv")

# in this case wa are going to use the individual samples saved from the script 20_object_WM_CX_processing
folder <- "../../data/individual_samples/"
sample_file <- dir(folder)

# set them up in a list
data.list <- lapply(sample_file, function(x){
  readRDS(paste0(folder,x))
}) %>%
  setNames(sample_file)

# save the list of features for the integration
features <- SelectIntegrationFeatures(object.list = data.list)
combined.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)

# this command creates an 'integrated' data assay
data.combined <- IntegrateData(anchorset = combined.anchors)

# score the cell cycle per 
DefaultAssay(data.combined) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data.combined <- CellCycleScoring(data.combined, s.features = s.genes, g2m.features = g2m.genes)
# head(data.combined@meta.data)
DefaultAssay(data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
# data.combined <- ScaleData(data.combined,vars.to.regress = c("percent.mt","nCount_RNA"))
data.combined <- ScaleData(data.combined,vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"))
# data.combined <- ScaleData(data.combined,vars.to.regress = c("percent.mt","nCount_RNA"))
# data.combined <- ScaleData(data.combined)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 0.2)

saveRDS(data.combined,"../../out/object/ManualClean/data.combined_WM_CX.rds")