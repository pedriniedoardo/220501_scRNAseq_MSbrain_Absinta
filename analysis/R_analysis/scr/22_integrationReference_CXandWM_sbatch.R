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

# pull the index of the samples
df_sample <- data.frame(file = sample_file) %>% 
  mutate(index = 1:nrow(.)) %>% 
  mutate(official_id = str_remove_all(file,pattern = ".rds"))

# seelct the samples that I want to use as reference
index_ref_sample <-  df_sample %>% 
  filter(official_id %in% c("s43","s57","s53","s39","s49","s60","s9","s3","s1","s27","s8","s23")) %>% 
  pull(index)

c("s43","s57","s53","s39","s49","s60","s9","s3","s1","s27","s8","s23")
df_sample %>% 
  filter(official_id %in% c("s43","s57","s53","s39","s49","s60","s9","s3","s1","s27","s8","s23"))

# set them up in a list
data.list <- lapply(sample_file, function(x){
  readRDS(paste0(folder,x))
}) %>%
  setNames(sample_file)

# make sure all the samples have been processed the same way: normalize and identify variable features for each dataset independently
data.list <- pmap(list(obj = data.list,name =names(data.list)), function(obj,name) {
  print(name)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  return(obj)
})

# save the list of features for the integration
features <- SelectIntegrationFeatures(object.list = data.list)

# make sure the PCA is available for each object
data.list <- pmap(list(obj = data.list,name =names(data.list)), function(obj,name) {
  print(name)
  obj <- ScaleData(obj, features = features, verbose = T) %>% 
    RunPCA(features = features, verbose = T)
  return(obj)
})

# use the reference approach with the defined samples
combined.anchors <- FindIntegrationAnchors(object.list = data.list,
                                           anchor.features = features,
                                           reference = index_ref_sample,
                                           reduction = "rpca")

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

saveRDS(data.combined,"../../out/object/ManualClean/data.combined_WM_CX_rpca.rds")