library(Seurat)
library(tidyverse)

data.combined <- readRDS("../../out/object/data.combined_fix_filter_norm_doublet_integrated_SoupX.rds")

DimPlot(data.combined,raster = T,split.by = "official_id",ncol=5)

data.combined@meta.data %>% dim()
data.combined@reductions$umap@cell.embeddings %>% dim()

test <- subset(data.combined,orig.ident=="11_69_control_cortex_martina_ref_CR6_1")

DimPlot(test)
