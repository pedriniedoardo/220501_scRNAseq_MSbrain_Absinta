# AIM ---------------------------------------------------------------------
# try to run harmony by merging the matrices from the individula objects. this will allow the skipping of the regular integration. the regular integration is needed as to run harmony the matrices should have the same number/order of the genes.

# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
library(harmony)
library(ggrepel)
library(ComplexHeatmap)

# read in the preliminary object ------------------------------------------
scobj_ref <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmony_test.rds")
scobj_ref$barcodes <- colnames(scobj_ref)

# select a random sample of 20k cells
set.seed(123)
sample_id <- sample(colnames(scobj_ref),size = 20000,replace = F)

# for the test try to downssample the dataset
# scobj_test <- subset(scobj_ref,downsample=500)
scobj_test <- subset(scobj_ref,barcodes %in% sample_id)

meta <- scobj_test@meta.data

# Run Harmony -------------------------------------------------------------
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.
# sobj_total_h <- scobj_test %>%
#   RunHarmony("orig.ident", plot_convergence = TRUE)

# sobj_total_h <- scobj_test %>%
#   RunHarmony("origin", plot_convergence = TRUE)

sobj_total_h <- scobj_test %>%
  RunHarmony(c("orig.ident","origin"), plot_convergence = TRUE)
# Harmony with two or more covariates
# Do the same with your Seurat object:
# seuratObject <- RunHarmony(seuratObject, c("dataset", "donor", "batch_id"))
# To directly access the new Harmony embeddings, use the Embeddings command.
harmony_embeddings <- Embeddings(sobj_total_h, 'harmony')
harmony_embeddings[1:5, 1:5]
# Let's make sure that the datasets are well integrated in the first 2 dimensions after Harmony.
# DimPlot(object = sobj_total_h, reduction = "harmony", pt.size = .1, group.by = "sample_id")

# Downstream analysis -----------------------------------------------------
# Many downstream analyses are performed on low dimensional embeddings, not gene expression. To use the corrected Harmony embeddings rather than PCs, set reduction = 'harmony'. For example, let's perform the UMAP and Nearest Neighbor analyses using the Harmony embeddings.
sobj_total_h <- sobj_total_h %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.2) %>%
  identity()

# verify that all the relevant slots are filled
sobj_total_h@assays$RNA@counts[1:20,1:10]
sobj_total_h@assays$RNA@data[1:20,1:10]
sobj_total_h@assays$RNA@scale.data[1:20,1:10]

dim(sobj_total_h@assays$RNA@counts)
dim(sobj_total_h@assays$RNA@data)
dim(sobj_total_h@assays$RNA@scale.data)

# DimPlot(sobj_total_h,group.by = "origin")
DimPlot(sobj_total_h,split.by = "origin")

DimPlot(sobj_total_h,split.by = "origin", group.by = "seurat_clusters_ref",label = T)

# saveRDS(sobj_total_h,"../out_large/scRNAseq_analysis/object/sobj_total_h_fix_filter_norm_doublet_harmony_5K.rds")
# saveRDS(sobj_total_h,"../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegration.rds")
