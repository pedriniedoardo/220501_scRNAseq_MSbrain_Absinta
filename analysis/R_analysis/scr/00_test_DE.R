# # libraries ---------------------------------------------------------------
# library(Seurat)
# library(SeuratData)
# library(tidyverse)
# library(future)
# library(SeuratWrappers)
# 
# # read in the data --------------------------------------------------------
# # read in the full object
# so <- readRDS("../../out/object/data.combined_fix_filter_norm_doublet_integrated_SoupX_harmony.rds")
# 
# # set the seed
# set.seed(1234)
# # define the number of cells to keep
# cell_keep <- 40000
# # define the index of the cells to keep
# id <- sample(c(rep(1,cell_keep),rep(0,(dim(so@meta.data)[1]-cell_keep))))
# # add the index to the metadata
# so$test <- id
# 
# # subset the dataset to the specified size
# test <-subset(so,subset = test==1)
# 
# # regular run -------------------------------------------------------------
# # data
# DefaultAssay(test) <- "RNA"
# Idents(test)<-"seurat_clusters"
# 
# # find markers for every cluster compared to all remaining cells, report only the positive
# # ones
# t <- Sys.time()
# out1 <- FindMarkers(object = test,ident.1 = "1",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Sys.time()-t
# 
# # use presto --------------------------------------------------------------
# t <- Sys.time()
# out2 <- RunPresto(object = test,ident.1 = "1",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Sys.time()-t
# 
# # use future plan ---------------------------------------------------------
# # library(future)
# # plan before planning
# plan()
# # change the planning
# plan("multisession", workers = 16)
# # confirm the planning
# plan()
# 
# t <- Sys.time()
# out3 <- FindMarkers(object = test,ident.1 = "1",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Sys.time()-t
# 
