# AIM ---------------------------------------------------------------------
# perform the new cleaning after we decided to remove the WM specific strange OLIGO cluster.

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
# read in the original dataset from the last version round of the analysis
sobj_ref <- readRDS("../../out/object/revision/112_WMCX_ManualClean2_harmonySkipIntegration_AllSoupX_4000.rds")

# filtering barcodes ------------------------------------------------------
# confirm the identity of the barcodes to be removed
DimPlot(sobj_ref,raster = T,label = T,group.by = "RNA_snn_res.0.3")

# Martina decided to remove all the barcodes associated to the cluster 7 and 11 from the resolution 0.3
whitelist_ref <- sobj_ref@meta.data %>%
  filter(RNA_snn_res.0.3 %in% c("7","11")) %>%
  rownames_to_column()
# number of cells to remove
dim(whitelist_ref)
# save the table
whitelist_ref %>%
  write_tsv("../../out/table/revision/114_whitelist_ref.tsv")

# in this case all the cells are already present in the sobj_ref, therefore I can filter from that one
# total cells
dim(sobj_ref)[2]

# cells to be removed
dim(whitelist_ref)[1]

# cells to retain after filtering
dim(sobj_ref)[2] - (dim(whitelist_ref)[1])

# perform the filtering
# add the metadata for the filteirng
meta_ref <- sobj_ref@meta.data %>%
  rownames_to_column()

# pull the barcodes from the whitelist and add the whitelist column to the metadata
# length(whitelist_NEU$rowname)
# sum(meta_ref$rowname %in% whitelist_NEU$rowname)

# length(whitelist_ref$rowname)
# sum(meta_ref$rowname %in% whitelist_ref$rowname)

# define the filtering variable
meta_ref_full <- meta_ref %>%
  mutate(whitelist = case_when(rowname %in% c(whitelist_ref$rowname)~0,
                               T~1))

# add the withelist ot the meta
sobj_ref$whitelist <- meta_ref_full$whitelist

# sanity check
sum(sobj_ref$whitelist)
dim(sobj_ref)[2] - (dim(whitelist_ref)[1])

# perform the filtering
sobj_filter <- subset(sobj_ref,subset = whitelist == 1)

# sanity check
dim(sobj_filter)[2]
dim(sobj_ref)[2] - (dim(whitelist_ref)[1])

DimPlot(sobj_filter,raster = T,label = T,group.by = "RNA_snn_res.0.3",split.by="origin")

# save the cleaned object -------------------------------------------------
saveRDS(sobj_filter,"../../out/object/revision/114_sobj_filter.rds")

