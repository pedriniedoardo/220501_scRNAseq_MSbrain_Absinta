# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(tidyverse)

library(cacoa)
library(conos)
library(sccore)

# sample process the seurat object ----------------------------------------
# # run the analysis to produce a sample PBMCs dataset
# # install dataset
# # InstallData("ifnb")
# # load dataset
# # LoadData("ifnb")
# load("data/PBMC/ifnb.SeuratData/data/ifnb.rda")
# # DimPlot(ifnb)
# # split the dataset into a list of two seurat objects (stim and CTRL)
# ifnb.list <- SplitObject(ifnb, split.by = "stim")
#
# # normalize and identify variable features for each dataset independently
# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
#
# # select features that are repeatedly variable across datasets for integration
# features <- SelectIntegrationFeatures(object.list = ifnb.list)
#
# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# # this command creates an 'integrated' data assay
# immune.combined <- IntegrateData(anchorset = immune.anchors)
#
# # Perform an integrated analysis
# # Now we can run a single integrated analysis on all cells!
# # specify that we will perform downstream analysis on the corrected data note that the
# # original unmodified data still resides in the 'RNA' assay
# DefaultAssay(immune.combined) <- "integrated"
#
# # Run the standard workflow for visualization and clustering
# immune.combined <- ScaleData(immune.combined, verbose = FALSE)
# immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
# immune.combined <- FindClusters(immune.combined, resolution = 0.5)
# # Visualization
# p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
# p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
# p1 + p2
# # scale also the RNA slot of the obeject
# DefaultAssay(immune.combined) <- "RNA"
# # Run the standard workflow for visualization and clustering
# immune.combined <- ScaleData(immune.combined, verbose = FALSE)
#
# # save the object as peocessed
# saveRDS(immune.combined,"out/object/seurat_PBMC.rds")

# build the cacoa object from seurat --------------------------------------
# load the seurat object

# all the cells
# so <- readRDS("/home/edo/Desktop/share/BS/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
so <- readRDS("../out_large/scRNAseq_analysis/object/sobj_total_h_fix_filter_norm_doublet_harmony_5K_scaleSbatch_SCtypeAnnotation.rds")
DimPlot(so,label = T,group.by = "annotation_confident")

# load the annotation for the subclusters
meta_subNeu <- read_tsv("out/table/meta_subNeu_refHippo_classification.tsv") %>%
  select(barcodes,contains("assigned"))

# wrangling ---------------------------------------------------------------
# add the annotation to the dataset
meta_ref <- so@meta.data %>%
  rownames_to_column("barcodes") %>%
  left_join(meta_subNeu,by = "barcodes") %>%
  mutate(assigned.predicted.id.l1 = case_when(is.na(assigned.predicted.id.l1)~annotation_confident,
                                              T~assigned.predicted.id.l1),
         assigned.predicted.id.l2 = case_when(is.na(assigned.predicted.id.l2)~annotation_confident,
                                              T~assigned.predicted.id.l2),
         assigned.predicted.id.l3 = case_when(is.na(assigned.predicted.id.l3)~annotation_confident,
                                              T~assigned.predicted.id.l3)) %>%
  column_to_rownames("barcodes")

# update the meta
so@meta.data <- meta_ref

# save the new heatmap
DimPlot(so,group.by = "assigned.predicted.id.l1",label = T)


# subset the dataset ------------------------------------------------------
so_subset <- subset(so,ADNPscoring%in%c("control","MCI-lowAD"))
so_subset <- subset(so,ADNPscoring%in%c("control","MCI-lowAD")&orig.ident.cca != "10_cr_61")
so_subset <- subset(so,ADNPscoring%in%c("control","MCI-lowAD")&orig.ident.cca != "03_cr_61")
table(so_subset$orig.ident.cca)

# create a fake sample varibale to run the tool. generally there should be replicates in the dataset
meta <- so_subset@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(rowname = cell) %>%
  column_to_rownames("rowname")

# put a sample code to build the con object from a seurat object
# Cacoa currently supports inputs in several formats (see below). Most of them require the following metadata:
#
# sample.groups: vector with condition labels per sample named with sample ids
sample.groups <- meta$ADNPscoring
names(sample.groups) <- meta$orig.ident.cca

# cell.groups: cell type annotation vector named by cell ids
cell.groups <- meta$assigned.predicted.id.l1
names(cell.groups) <- meta$cell

# sample.per.cell: vector with sample labels per cell named with cell ids
sample.per.cell <- meta$orig.ident.cca
names(sample.per.cell) <- meta$cell

# ref.level: id of the condition, corresponding to the reference (i.e. control)
table(sample.groups)
ref.level <- "control"

# target.level: id of the condition, corresponding to the target (i.e. case)
target.level <- "MCI-lowAD"

# Additionally, embedding parameter containing a matrix or data.frame with a cell embedding can be provided. Rownames should match to the cell ids. It is used for visualization and some cluster-free analysis.
embedding <- so_subset@reductions$umap@cell.embeddings %>%
  data.frame()

# Parameter graph.name is required for cluster-free analysis, and must contain a name of joint graph in Seurat object. For that, the Seurat object must have a joint graph estimated (see FindNeighbors). For visualization purposes, Seurat also must have cell embedding estimated or the embedding data frame must be provided in the embedding parameter.
graph.name <- "integrated_nn"

# Seurat object so
cao <- Cacoa$new(so_subset,
                 sample.groups=sample.groups,
                 cell.groups=cell.groups,
                 sample.per.cell=sample.per.cell,
                 ref.level=ref.level,
                 target.level=target.level,
                 embedding = embedding,
                 graph.name=graph.name)

# run cacoa ---------------------------------------------------------------

# Cluster-based changes ---------------------------------------------------
# The fastest way to visualize changes in the dataset is to show, what cell types changed their abundance and expression patterns.

# Estimate cluster-based changes
cao$estimateCellLoadings()
cao$estimateExpressionShiftMagnitudes()

# Plot compositional changes
cao$plotCellLoadings(show.pvals=FALSE)
# save plot in case all the sample are used
# ggsave("out/image/cacoa_MCIvsCTRL.pdf",width = 10,height = 6)

# The red line here shows statistical significance.
# Plot expression changes
cao$plotExpressionShiftMagnitudes()

# Here, y-axis shows magnitude of changes, cells show both expression and composition changes.


# subset control vs NOLD --------------------------------------------------
so_subset2 <- subset(so,ADNPscoring%in%c("control","NOLD-lowAD"))
so_subset2 <- subset(so,ADNPscoring%in%c("control","NOLD-lowAD")&orig.ident.cca != "10_cr_61")
so_subset2 <- subset(so,ADNPscoring%in%c("control","NOLD-lowAD")&orig.ident.cca != "03_cr_61")
table(so_subset2$orig.ident.cca)

# create a fake sample varibale to run the tool. generally there should be replicates in the dataset
meta2 <- so_subset2@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(rowname = cell) %>%
  column_to_rownames("rowname")

# put a sample code to build the con object from a seurat object
# Cacoa currently supports inputs in several formats (see below). Most of them require the following metadata:
#
# sample.groups: vector with condition labels per sample named with sample ids
sample.groups2 <- meta2$ADNPscoring
names(sample.groups2) <- meta2$orig.ident.cca

# cell.groups: cell type annotation vector named by cell ids
cell.groups2 <- meta2$assigned.predicted.id.l1
names(cell.groups2) <- meta2$cell

# sample.per.cell: vector with sample labels per cell named with cell ids
sample.per.cell2 <- meta2$orig.ident.cca
names(sample.per.cell2) <- meta2$cell

# ref.level: id of the condition, corresponding to the reference (i.e. control)
table(sample.groups2)
ref.level2 <- "control"

# target.level: id of the condition, corresponding to the target (i.e. case)
target.level2 <- "NOLD-lowAD"

# Additionally, embedding parameter containing a matrix or data.frame with a cell embedding can be provided. Rownames should match to the cell ids. It is used for visualization and some cluster-free analysis.
embedding2 <- so_subset2@reductions$umap@cell.embeddings %>%
  data.frame()

# Parameter graph.name is required for cluster-free analysis, and must contain a name of joint graph in Seurat object. For that, the Seurat object must have a joint graph estimated (see FindNeighbors). For visualization purposes, Seurat also must have cell embedding estimated or the embedding data frame must be provided in the embedding parameter.
graph.name <- "integrated_nn"

# Seurat object so
cao2 <- Cacoa$new(so_subset2,
                  sample.groups=sample.groups2,
                  cell.groups=cell.groups2,
                  sample.per.cell=sample.per.cell2,
                  ref.level=ref.level2,
                  target.level=target.level2,
                  embedding = embedding2,
                  graph.name=graph.name)

# run cacoa ---------------------------------------------------------------

# Cluster-based changes ---------------------------------------------------
# The fastest way to visualize changes in the dataset is to show, what cell types changed their abundance and expression patterns.

# Estimate cluster-based changes
cao2$estimateCellLoadings()
cao2$estimateExpressionShiftMagnitudes()

# Plot compositional changes
cao2$plotCellLoadings(show.pvals=FALSE)
# save plot in case all the sample are used
# ggsave("out/image/cacoa_NOLDvsCTRL.pdf",width = 10,height = 6)

# The red line here shows statistical significance.
# Plot expression changes
cao2$plotExpressionShiftMagnitudes()

