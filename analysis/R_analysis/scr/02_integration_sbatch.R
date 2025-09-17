# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)

# read in the data --------------------------------------------------------
# load the LUT
LUT <- read_csv("../../data/LUT_sample.csv")

# in this case wa are going to use the fix threshold filtered data
df_file <- dir("../../out/object/") %>%
  str_subset(pattern = c("datasc_fix_filter_norm_doublet_.*_SoupX")) %>%
  str_subset(pattern = c("Singlet"),negate = F) %>%
  data.frame(sample_file = .) %>%
  mutate(sample_id = str_remove_all(sample_file,pattern = "datasc_fix_filter_norm_doublet_Singlet_|_SoupX|.rds")) %>%
  left_join(.,LUT,by = c("sample_id"))

# set them up in a list
data.list <- lapply(df_file$sample_file, function(x){
  readRDS(paste0("../../out/object/",x))
}) %>%
  setNames(df_file$sample_id)

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

saveRDS(data.combined,"../../out/object/data.combined_fix_filter_norm_doublet_integrated_SoupX.rds")