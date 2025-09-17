# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(tidyverse)
library(future)

# read in the data --------------------------------------------------------
# above I did the generation of the object
so <- readRDS("../../out/object/ManualClean/data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony.rds")

# Identify conserved cell type markers ------------------------------------
# data
DefaultAssay(so) <- "RNA"
so@meta.data

Idents(so)<-"seurat_clusters"
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
t <- Sys.time()
# library(future)
# plan before planning
# plan()
# # change the planning
# plan("multisession", workers = 16)
# # confirm the planning
# plan()
# sobj_total_h.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sobj_total_h.markers <- RunPrestoAll(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Sys.time()-t

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/ManualClean/FindAllMarkers_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_PosOnly.tsv")

# save the top 100 markers per cluster
sobj_total_h.markers %>%
  group_by(cluster) %>%
  mutate(rank = rank(order(p_val_adj, -abs(avg_log2FC)), ties.method='first')) %>%
  arrange(cluster,rank) %>%
  filter(rank < 101) %>%
  write_tsv("../../out/table/ManualClean/FindAllMarkers_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_PosOnly_top100.tsv")

# save the average expression of each gene per cluster
DefaultAssay(so)
Idents(so)<-"seurat_clusters"

df_avg_exp <- AverageExpression(so)$RNA %>% 
  data.frame() %>% 
  rownames_to_column("genes")

dim(df_avg_exp)

# save the table of the average expressions
df_avg_exp %>% 
  write_tsv("../../out/table/ManualClean/AvgExp_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_SeuratClusters.tsv")
