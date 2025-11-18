# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the final object ------------------------------------------------
scobj <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")

DimPlot(scobj,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(scobj,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(scobj,label = T,raster = T,group.by = "expertAnno.l1",split.by = "origin")
DimPlot(scobj,label = T,raster = T,group.by = "pathology_class")

# wrangling ---------------------------------------------------------------
# explore the full meta
scobj@meta.data

# subset the dataset only for the new samples before being deposited
sobj_new <- subset(scobj,subset = orig.ident %in% c("s25","s26","s27","s31","s35","s39","s41","s42","s43","s44","s45","s46","s47","s48","s49","s51","s53","s54","s55","s57","s60","s65","s69"))

dim(sobj_new)
# export the metadata for the new samples
df_meta <- sobj_new@meta.data %>%
  rownames_to_column("barcodes")

df_meta %>%
  write_tsv("../../out/table/revision/metadata_GEOsubmission.tsv")

df_meta %>%
  group_by(orig.ident) %>%
  summarise(n = n()) %>%
  pull(n) %>%
  sum()

scobj@meta.data %>%
  group_by(orig.ident) %>%
  summarise(n = n()) %>% 
  # pull(orig.ident) %>%
  dplyr::filter(orig.ident %in% c("s25","s26","s27","s31","s35","s39","s41","s42","s43","s44","s45","s46","s47","s48","s49","s51","s53","s54","s55","s57","s60","s65","s69")) %>%
  pull(n) %>%
  sum()

# export the normalized reads
GetAssayData(sobj_new,assay = "RNA",slot = "data") %>%
  saveRDS("../../out/object/revision/normCount_GEOsubmission.rds")

# export the raw reads
GetAssayData(sobj_new,assay = "RNA",slot = "count") %>%
  saveRDS("../../out/object/revision/rawCount_GEOsubmission.rds")
