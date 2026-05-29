# AIM ---------------------------------------------------------------------
# custom analysis for correlation of TSPO with signature scores provided by Sofia

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

# read in the dataset -----------------------------------------------------
# read in the reference object
data.combined <- readRDS("../../out/object/129_MG_subcluster_HarmonySample_martinaCluster.rds")
DimPlot(data.combined,label = T,raster = F,group.by = "cluster_martina")

# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("TSPO")

table(data.combined@meta.data$cluster_martina)

# read in the signatures to analyze
