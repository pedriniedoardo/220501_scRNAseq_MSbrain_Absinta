# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(enrichR)
library(patchwork)


# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegAllSoupX_expertAnno.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")

# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("TSPO")

table(data.combined@meta.data$annotation_confident)

# add the SIT annotatation
# data.combined@meta.data$cell_type2 <- data.combined@meta.data |> 
#   mutate(cell_type2 = case_when(# according to label transfer cluster 12 and 6 are Bipolar cells
#     cell_type %in% c("BP_0","BP_10","BP_11","BP_13","BP_16","BP_17","BP_5","BP_8","12","6")~"BP",
#     T~cell_type)) |> 
#   pull(cell_type2)


# plot --------------------------------------------------------------------
# add the same grouping and style of the senmayo for the summary of the SIT score
