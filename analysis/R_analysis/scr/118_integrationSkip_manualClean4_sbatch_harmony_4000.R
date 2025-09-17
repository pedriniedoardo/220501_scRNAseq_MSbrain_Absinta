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
# read in the cleaned dataset
sobj_filter <- readRDS("../../out/object/revision/117_sobj_filter.rds")

DimPlot(sobj_filter,raster = T,group.by = "seurat_clusters",label = T,split.by = "origin")

# wrangling ---------------------------------------------------------------
# save the new meta with only the info I want to track from the previous analysis
# meta_sobj_filter <- sobj_filter@meta.data %>%
#   rownames_to_column("barcodes") %>%
#   select("barcodes",
#          "orig.ident",
#          "seurat_clusters_ref",
#          "sample_id",
#          "facility",
#          "origin",
#          "disease",
#          "sex",
#          "age",
#          "pathology",
#          "pathology_class",
#          "patient",
#          "plaque",
#          "PMI",
#          "seurat_clusters.old",
#          "origin_alt",
#          "annotation_confident.old",
#          "expertAnno.l1.old",
#          "expertAnno.l2.old") %>%
#   column_to_rownames("barcodes")

meta_sobj_filter <- sobj_filter@meta.data %>%
  rownames_to_column("barcodes") %>%
  select("barcodes",
         "orig.ident",
         "seurat_clusters_ref",
         "sample_id",
         "facility",
         "origin",
         "disease",
         "sex",
         "age",
         "pathology",
         "pathology_class",
         "patient",
         "plaque",
         "PMI",
         seurat_clusters.old = "seurat_clusters",
         "origin_alt",
         annotation_confident.old = "annotation_confident",
         expertAnno.l1.old = "expertAnno.l1",
         expertAnno.l2.old = "expertAnno.l2") %>%
  column_to_rownames("barcodes")


# sc processing -----------------------------------------------------------
# I need to create a single object to add the cell cycle scoring and other metadata. in this case do not perform any new filtering
sobj_total <- CreateSeuratObject(counts = sobj_filter@assays$RNA@counts,
                                 meta.data = meta_sobj_filter,
                                 project = "CXWM_ManualClean4",
                                 min.cells = 0,
                                 min.features = 0)

# clean up the memory
remove(list = str_subset(ls(),pattern = "sobj_total",negate = T))
gc()

# add the cell cycle analysis
DefaultAssay(sobj_total) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes, g2m.features = g2m.genes)
sobj_total$percent.mt <- PercentageFeatureSet(sobj_total, pattern = "^MT-")
sobj_total$percent.ribo <- PercentageFeatureSet(sobj_total, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
sobj_total$percent.globin <- Seurat::PercentageFeatureSet(sobj_total,pattern = "^HB[^(P)]")

# pull all the variable to scale
sobj_total@assays$RNA@scale.data
# all.genes <- rownames(sobj_total)

# rescale the data for regressing out the sources of variation do not scale all the genes. if needed scale them before the heatmap call
sobj_total <- sobj_total %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>%
  # I can scale the missing features afterwards now focus on the highly variable one for speed purposes
  ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score","origin","facility"), verbose = T) %>% 
  # ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T,features = all.genes) %>% 
  RunPCA(npcs = 30, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:30,return.model = TRUE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  identity()

# Run Harmony -------------------------------------------------------------
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.
sobj_total_h <- sobj_total %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
# Harmony with two or more covariates
# Do the same with your Seurat object:
# seuratObject <- RunHarmony(seuratObject, c("dataset", "donor", "batch_id"))
# To directly access the new Harmony embeddings, use the Embeddings command.
harmony_embeddings <- Embeddings(sobj_total_h, 'harmony')
harmony_embeddings[1:5, 1:5]
# Let's make sure that the datasets are well integrated in the first 2 dimensions after Harmony.
# DimPlot(object = sobj_total_h, reduction = "harmony", pt.size = .1, group.by = "sample_id")

# Downstream analysis -----------------------------------------------------
# Many downstream analyses are performed on low dimensional embeddings, not gene expression. To use the corrected Harmony embeddings rather than PCs, set reduction = 'harmony'. For example, let's perform the UMAP and Nearest Neighbor analyses using the Harmony embeddings.
sobj_total_h <- sobj_total_h %>%
  RunUMAP(reduction = "harmony", dims = 1:30,return.model = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  identity()

# verify that all the relevant slots are filled
sobj_total_h@assays$RNA@counts[1:20,1:10]
sobj_total_h@assays$RNA@data[1:20,1:10]
sobj_total_h@assays$RNA@scale.data[1:20,1:10]

dim(sobj_total_h@assays$RNA@counts)
dim(sobj_total_h@assays$RNA@data)
dim(sobj_total_h@assays$RNA@scale.data)

# DimPlot(sobj_total_h,group.by = "origin")
# DimPlot(sobj_total_h,split.by = "origin")

# save the final object ---------------------------------------------------
# saveRDS(sobj_total_h,"../out_large/scRNAseq_analysis/object/sobj_total_h_fix_filter_norm_doublet_harmony_5K.rds")
saveRDS(sobj_total_h,"../../out/object/revision/118_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000.rds")
