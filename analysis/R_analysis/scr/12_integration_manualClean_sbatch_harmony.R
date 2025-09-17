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
# read in the integrated object
combined <- readRDS("../../out/object/ManualClean/data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean.rds")

# save the meta of the combined object
df_meta <- combined@meta.data %>%
  rownames_to_column("barcode")

# # combined@reductions$umap@cell.embeddings
# 
# plot <- DimPlot(object = combined,group.by = "seurat_clusters")
# plot$data %>% 
#   ggplot(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters))+geom_point(size=0.3)

# split the object in independent datasets --------------------------------
# split the integrated objects
list_split <- SplitObject(combined, split.by = "orig.ident")
# remove the original object to make the session slimmer
remove(combined)
gc()

# save the sperse matrices for each object
total_sm <- lapply(list_split,function(x){
  x@assays$RNA@counts
}) %>%
  do.call(what = cbind,args = .)

# check the total_sm
dim(total_sm)
# confirmt he total number of cells

# lapply(list_split,function(x){
#   dim(x@assays$RNA@counts)
# })

lapply(list_split,function(x){
  dim(x@assays$RNA@counts)[2]
}) %>%
  unlist() %>%
  sum()
# check a sample fo the dataset
total_sm[1:10,50:70]

# notice it is critacal that both matrices have the same dimension
# create a unique seurat object
# sobj_total <- CreateSeuratObject(counts = total_sm, project = "snAD", min.cells = 20, min.features = 200) %>%
#   Seurat::NormalizeData(verbose = T) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = T) %>%
#   RunPCA(pc.genes = pbmc@var.genes, npcs = 30, verbose = T) %>%
#   RunUMAP(reduction = "pca", dims = 1:30) %>%
#   FindNeighbors(reduction = "pca", dims = 1:30) %>%
#   FindClusters(resolution = 0.2) %>%
#   identity()
# I need to create a single object to add the cell cycle scoring and other metadata
sobj_total <- CreateSeuratObject(counts = total_sm, project = "snCx", min.cells = 20, min.features = 200)

# add the cell cycle analysis
DefaultAssay(sobj_total) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes, g2m.features = g2m.genes)
sobj_total$percent.mt <- PercentageFeatureSet(sobj_total, pattern = "^MT-")
sobj_total$percent.ribo <- PercentageFeatureSet(sobj_total, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

# add all the original metadata 
meta_new <- sobj_total@meta.data %>% 
  rownames_to_column("barcode")

meta_new_total <- left_join(meta_new,df_meta,"barcode",suffix=c(".harmony",".cca")) %>% 
  column_to_rownames("barcode")

# update the meta
sobj_total@meta.data <- meta_new_total

# pull all the variable to scale
sobj_total@assays$RNA@scale.data
# all.genes <- rownames(sobj_total)

# rescale the data for regressing out the sources of variation
sobj_total <- sobj_total %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  # I can scale the missing features afterwards
  ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T) %>% 
  # ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T,features = all.genes) %>% 
  RunPCA(npcs = 30, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2) %>%
  identity()

# Run Harmony -------------------------------------------------------------
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.
sobj_total_h <- sobj_total %>%
  RunHarmony("orig.ident.cca", plot_convergence = TRUE)
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
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.2) %>%
  identity()

# verify that all the relevant slots are filled
sobj_total_h@assays$RNA@counts[1:20,1:10]
sobj_total_h@assays$RNA@data[1:20,1:10]
sobj_total_h@assays$RNA@scale.data[1:20,1:10]

dim(sobj_total_h@assays$RNA@counts)
dim(sobj_total_h@assays$RNA@data)
dim(sobj_total_h@assays$RNA@scale.data)

# saveRDS(sobj_total_h,"../out_large/scRNAseq_analysis/object/sobj_total_h_fix_filter_norm_doublet_harmony_5K.rds")
saveRDS(sobj_total_h,"../../out/object/ManualClean/data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony.rds")
