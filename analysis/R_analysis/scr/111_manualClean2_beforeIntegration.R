# AIM ---------------------------------------------------------------------
# perform the new cleaning after the exploration of the NEU subcluster

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
sobj_ref <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegAllSoupX_expertAnno.rds")

# read in the subcluster NEU data
sobj_NEU <- readRDS("../../out/object/revision/100_subcluster_NEU.rds")

# read in the subcluster OLIGO data
sobj_OLIGO <- readRDS("../../out/object/revision/100_subcluster_OLIGO.rds")

# filtering barcodes ------------------------------------------------------
# confirm the identity of the barcodes to be removed
DimPlot(sobj_NEU,raster = T,label = T,group.by = "RNA_snn_res.0.2")
# Martina decided to remove all the barcodes associated to the cluster 0 from RNA_snn_res.0.2
whitelist_NEU <- sobj_NEU@meta.data %>%
  filter(RNA_snn_res.0.2 %in% c("0","14")) %>%
  rownames_to_column()
# number of cells to remove
dim(whitelist_NEU)
# save the table
whitelist_NEU %>%
  write_tsv("../../out/table/revision/111_whitelist_NEU.tsv")

# confirm the identity of the barcodes to be removed
DimPlot(sobj_OLIGO,raster = T,label = T,group.by = "RNA_snn_res.0.2")
# Martina decided to remove all the barcodes associated to the cluster 0 from RNA_snn_res.0.2
whitelist_OLIGO <- sobj_OLIGO@meta.data %>%
  filter(RNA_snn_res.0.2 %in% c("4","7")) %>%
  rownames_to_column()
# number of cells to remove
dim(whitelist_OLIGO)
# save the table
whitelist_OLIGO %>%
  write_tsv("../../out/table/revision/111_whitelist_OLIGO.tsv")

# clean also from the original sample
DimPlot(sobj_ref,raster = T,label = T,group.by = "seurat_clusters")
# Martina decided to remove all the barcodes associated to cluster 14 from seurat_clusters metadata
whitelist_ref <- sobj_ref@meta.data %>%
  filter(seurat_clusters %in% c("14")) %>%
  rownames_to_column()
# number of cells to remove
dim(whitelist_ref)
# save the table
whitelist_ref %>%
  write_tsv("../../out/table/revision/111_whitelist_ref.tsv")

# in this case all the cells are already present in the sobj_ref, therefore I can filter from that one
# total cells
dim(sobj_ref)[2]

# cells to be removed
dim(whitelist_ref)[1]+dim(whitelist_NEU)[1]+dim(whitelist_OLIGO)[1]

# cells to retain after filtering
dim(sobj_ref)[2] - (dim(whitelist_ref)[1]+dim(whitelist_NEU)[1]+dim(whitelist_OLIGO)[1])

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
  mutate(whitelist = case_when(rowname %in% c(whitelist_OLIGO$rowname,whitelist_NEU$rowname,whitelist_ref$rowname)~0,
                               T~1))

# add the withelist ot the meta
sobj_ref$whitelist <- meta_ref_full$whitelist

# sanity check
sum(sobj_ref$whitelist)
dim(sobj_ref)[2] - (dim(whitelist_ref)[1]+dim(whitelist_NEU)[1]+dim(whitelist_OLIGO)[1])

# perform the filtering
sobj_filter <- subset(sobj_ref,subset = whitelist == 1)

# sanity check
dim(sobj_filter)[2]
dim(sobj_ref)[2] - (dim(whitelist_ref)[1]+dim(whitelist_NEU)[1]+dim(whitelist_OLIGO)[1])

# save the cleaned object -------------------------------------------------
saveRDS(sobj_filter,"../../out/object/revision/111_sobj_filter.rds")

