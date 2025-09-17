# AIM ---------------------------------------------------------------------
# run the subcluster analysis for the AST cells from the WM+CX dataset
# this one will specifically filter out the cells that martina wanted to remove and annotate using her annotation

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

# read in the data --------------------------------------------------------
# read in the cells from run 02
sobj <- readRDS("../../out/object/131_AST_subcluster_HarmonySample.rds")
DimPlot(sobj,raster = T,group.by = "RNA_snn_res.0.5",label = T)
DimPlot(sobj,raster = T,group.by = "RNA_snn_res.0.5",label = T,split.by = "origin")

# wrangling ---------------------------------------------------------------
# subset only the cells of interest.
# martina asked not to rerun the analysis after removeing the cells, but just to remove them from the current object
id_keep <- data.frame(id = levels(sobj$RNA_snn_res.0.5)) %>%
  filter(!(id %in% c("9","10","11","12") )) %>%
  pull(id)

sobj_subset <- subset(sobj,subset = RNA_snn_res.0.5 %in% id_keep)
DimPlot(sobj_subset,group.by = "RNA_snn_res.0.5",label=T,raster = T)

# add the new annotatition she recommended
df_meta <- sobj_subset@meta.data

sobj_subset$cellid <- df_meta %>%
  mutate(cellid = case_when(RNA_snn_res.0.5 %in% c(1,2)~"Cortical AST",
                            RNA_snn_res.0.5 %in% c(3,4,7)~"Stressed",
                            RNA_snn_res.0.5 %in% c(8,14)~"Reactive",
                            RNA_snn_res.0.5 %in% c(0)~"Homeo",
                            RNA_snn_res.0.5 %in% c(13)~"AIMS",
                            RNA_snn_res.0.5 %in% c(6)~"Ciliated",
                            RNA_snn_res.0.5 %in% c(5)~"Perinodal")) %>%
  pull(cellid)

DimPlot(sobj_subset,raster = T,group.by = "cellid",label = T)

# sobj_total_h <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_01000_6000_15.rds")
saveRDS(sobj_subset,"../../out/object/132_AST_subcluster_HarmonySample_martina.rds")
