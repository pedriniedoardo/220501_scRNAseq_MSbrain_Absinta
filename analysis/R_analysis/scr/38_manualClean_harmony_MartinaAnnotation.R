# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(SeuratWrappers)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegAllSoupX_SCtypeAnnotation.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "annotation_confident")

# add the costume annotation defined by Martina
meta <- data.combined@meta.data

# create a LUT for the annotation of the clusters
LUT_meta <- data.frame(seurat_clusters = 0:14) |> 
  mutate(expertAnno.l1 = case_when(seurat_clusters %in% c(5)~"MG",
                                   seurat_clusters %in% c(13)~"LYM",
                                   seurat_clusters %in% c(0,4,14)~"OLIG",
                                   seurat_clusters %in% c(9)~"OPC",
                                   seurat_clusters %in% c(3,12)~"ASTRO",
                                   seurat_clusters %in% c(11)~"VAS",
                                   seurat_clusters %in% c(1)~"NEU_PY",
                                   seurat_clusters %in% c(2,6,10)~"NEU_EXC",
                                   seurat_clusters %in% c(7,8)~"NEU_INH"),
         expertAnno.l2 = case_when(seurat_clusters %in% c(5)~"MG",
                                   seurat_clusters %in% c(13)~"LYM",
                                   seurat_clusters %in% c(0,4,14)~"OLIG",
                                   seurat_clusters %in% c(9)~"OPC",
                                   seurat_clusters %in% c(3,12)~"ASTRO",
                                   seurat_clusters %in% c(11)~"VAS",
                                   seurat_clusters %in% c(1)~"NEU_PY",
                                   seurat_clusters %in% c(2)~"NEU_EXC_UPPER",
                                   seurat_clusters %in% c(6,10)~"NEU_EXC_LOWER",
                                   seurat_clusters %in% c(7,8)~"NEU_INH")) |> 
  mutate(seurat_clusters = as.factor(seurat_clusters))

# add the new annotation to the full metadata
meta_full <- left_join(meta,LUT_meta,by = "seurat_clusters")

# update the object with the new annotation
data.combined$expertAnno.l1 <- meta_full$expertAnno.l1
data.combined$expertAnno.l2 <- meta_full$expertAnno.l2

DimPlot(data.combined,group.by = "expertAnno.l1",raster = T,label = T)
DimPlot(data.combined,group.by = "expertAnno.l2",raster = T,label = T)

# save a summary table with average count per cell
summary_meta_new <- data.combined@meta.data %>%
  group_by(orig.ident,expertAnno.l1,pathology_class) %>%
  summarise(n = n(),
            mean_nCount = mean(nCount_RNA),
            mean_nFeature = mean(nFeature_RNA)) %>%
  ungroup()%>%
  group_by(orig.ident) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

write_tsv(summary_meta_new,"../../out/table/summary_meta_expertAnno.tsv")

# save the final object
saveRDS(data.combined,"../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegAllSoupX_expertAnno.rds")


# run DE ------------------------------------------------------------------
data.combined <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegAllSoupX_expertAnno.rds")

# define the comparison of interest
Idents(data.combined) <- "expertAnno.l1"

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/FindAllMarkers_WM_CX_harmonySkipIntegAllSoupX_expertAnno.tsv")

# filter the top 100 genes per cell type
sobj_total_h.markers %>%
  group_by(cluster) %>%
  slice(1:100) %>%
  write_tsv("../../out/table/FindAllMarkers_WM_CX_harmonySkipIntegAllSoupX_expertAnno_top100.tsv")

