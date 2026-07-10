# AIM ---------------------------------------------------------------------
# priduce some custom plot for sofia

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(cowplot)
library(patchwork)
library(Nebulosa)
library(pals)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
data.combined2 <- readRDS("../../data/all20_immune.rds")


# plot general ------------------------------------------------------------
# custom plot main population
DimPlot(data.combined,label = T,raster = T,group.by = "RNA_snn_res.0.1")

# Determine the number of clusters you need colors for
# RColorBrewer::display.brewer.all()
num_clusters <- length(levels(data.combined$RNA_snn_res.0.1))

# Generate a palette (e.g., "Set1" or "Paired")
# my_colors <- brewer.pal(n = num_clusters, name = "Paired")
color_id <- tableau20(num_clusters)
# names(color_id) <- unique(data.combined$identity %>% unlist())

show_col(color_id)

# Apply to DimPlot
p01 <- DimPlot(data.combined, 
        label = TRUE, 
        raster = TRUE, 
        group.by = "RNA_snn_res.0.1", 
        cols = color_id)
ggsave(plot = p01,filename = "../../out/image/revision/00_UMAP_brain_general_paletteUpdate.pdf",width = 6,height = 5)

# plot Absinta2021 --------------------------------------------------------
# custom plot main population
DimPlot(data.combined2,label = T,raster = T,group.by = "seurat_clusters")

# Determine the number of clusters you need colors for
# RColorBrewer::display.brewer.all()
num_clusters2 <- length(levels(data.combined2$seurat_clusters))

# Generate a palette (e.g., "Set1" or "Paired")
# my_colors <- brewer.pal(n = num_clusters, name = "Paired")
color_id2 <- tableau20(num_clusters2)
# names(color_id) <- unique(data.combined$identity %>% unlist())

show_col(color_id2)

# Apply to DimPlot
p02 <- DimPlot(data.combined2, 
               label = F, 
               raster = F, 
               group.by = "seurat_clusters", 
               cols = color_id2)
ggsave(plot = p02,filename = "../../out/image/revision/00_UMAP_IMM_Absinta2021_paletteUpdate.pdf",width = 6,height = 5)

