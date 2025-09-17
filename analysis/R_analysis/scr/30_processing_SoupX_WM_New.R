# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
# library(DoubletFinder)
library(glmGamPoi)

# READ IN DATA ------------------------------------------------------------
# read in the new data
id_sample <- dir("../../data/SoupX_default_WM/") %>%
  str_subset(pattern = "s25|s26|s27|s31")

# test <- readRDS(paste0("../../data/WM_data/WM_data_Abisnta_new/",x,".rds"))

# read in the processed object and save the metadata
list_meta_WM_ref <-lapply(id_sample,function(x){
  obj <- readRDS(paste0("../../data/WM_data/WM_data_Abisnta_new/",x,".rds"))
  obj@meta.data %>%
    data.frame() %>%
    rownames_to_column("barcodes")
}) %>% 
  setNames(id_sample)

# order the samples in the ref metadata with the same order of the id_sample
list_meta_WM_ref_fix <- list_meta_WM_ref[id_sample] 

# in this case I am relying on the ref metadata, therefore I will skip all the calssical QC filtering and rely on the cells kept in the ref metadata
# x <- "s31"
list_datasc <- lapply(id_sample,function(x){
  print(x)
  data <- Read10X(data.dir = paste0("../../data/SoupX_default_WM/",x))
  
  # datasc <- CreateSeuratObject(counts = data, project = LUT %>%
  #                                filter(sample_id == x) %>%
  #                                pull(sample_id), min.cells = 20, min.features = 200)
  datasc <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
  
  # # datasc@meta.data
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  # datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  
  # # label the cells based on the mt reads content
  # datasc$mt_bin <- datasc@meta.data %>%
  #   mutate(test = case_when(percent.mt < 1~"low",
  #                           percent.mt < 10~"mid",
  #                           T ~ "high")) %>%
  #   pull(test)
  
  # all the other covariates in the LUT table
  # load the meta ref
  meta_ref <- list_meta_WM_ref_fix[[x]] %>% 
    # separate(barcodes,into = c("barcode_id","sample_num"),sep = "-",remove = F) %>% 
    dplyr::select(barcodes,seurat_clusters_ref=seurat_clusters,sample_id = orig.ident)
  
  # save the meta of this object and add the filtering variable
  meta_test <- datasc@meta.data %>% 
    rownames_to_column("barcodes") %>% 
    # separate(barcodes,into = c("barcode_id","sample_num"),sep = "-",remove = F) %>% 
    # dplyr::select(-sample_num) %>% 
    left_join(meta_ref,by = "barcodes") %>% 
    # add a filtering variable
    mutate(test = case_when(is.na(sample_id)~0,
                            T~1)) %>% 
    column_to_rownames("barcodes")
  
  # swap the metadata in the object
  datasc@meta.data <- meta_test
  
  # # add the filtering variable based on the fixed threshold
  # datasc$test <- datasc@meta.data %>%
  #   mutate(test = percent.mt < 10 & nFeature_RNA > 500 & nFeature_RNA < 7000) %>%
  #   pull(test)
  # 
  # # add the filtering variable based on the
  # stats <- cbind(log10(datasc@meta.data$nCount_RNA), log10(datasc@meta.data$nFeature_RNA),
  #                datasc@meta.data$percent.mt)
  # 
  # # library(robustbase)
  # outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
  # #library(scater)
  # multi.outlier <- isOutlier(outlying, type = "higher")
  # # summary(multi.outlier)
  # 
  # datasc$not_outlier <- !as.vector(multi.outlier)
  # 
  return(datasc)
}) %>%
  setNames(id_sample)

# plot QC -----------------------------------------------------------------
# this step is not modifine the object and is already run tin the 00_explore_fixed_threshold
# meta_total <- lapply(list_datasc, function(x){
#   x@meta.data
# }) %>%
#   bind_rows(.id = "dataset") %>%
#   rownames_to_column("barcode")

# filtering of the cells --------------------------------------------------
# perform the filtering based on the fixed threshold
list_datasc_fixed <- lapply(list_datasc,function(x){
  datasc_filter <- subset(x, subset = test == 1)
})

# conform that the size of the filtered objects and the size of the metadata are compatibles
lapply(list_meta_WM_ref_fix,function(x){
  dim(x)[1]
}) %>% 
  unlist()

lapply(list_datasc_fixed,function(x){
  dim(x)[2]
}) %>% 
  unlist()

# there is a samll difference that probably is due to the automatic filtering in cellranger. After chatting with Martina she decided to proceed anyway

# proprocess the dataset --------------------------------------------------
# use the Normalize function. for the doublet removla I need to run the full preprocessing steps
list_datasc_fixed_norm <- lapply(list_datasc_fixed, function(x){
  datasc_filter <- x
  datasc_filter <- NormalizeData(datasc_filter)
  datasc_filter <- FindVariableFeatures(datasc_filter, selection.method = "vst", nfeatures = 2000)
  datasc_filter <- ScaleData(datasc_filter) 
  datasc_filter <- RunPCA(datasc_filter) 
  datasc_filter <- FindNeighbors(datasc_filter,dims = 1:30) 
  datasc_filter <- FindClusters(datasc_filter) 
  datasc_filter <- RunUMAP(datasc_filter,dims = 1:30) 
  
  datasc_filter
})
# saveRDS(list_datasc_fixed_norm,"../../out/object/test_SCTransform/list_datasc_fixed_norm.rds")
# list_datasc_fixed_norm <- readRDS("../../out/object/test_SCTransform/list_datasc_fixed_norm.rds")

list_datasc_fixed_norm[[1]]@meta.data %>% colnames()

# subset the datasets after removal of the doublets befreo the integration
pmap(list(names(list_datasc_fixed_norm),list_datasc_fixed_norm),function(x,y){
  # fix the meta
  # meta <- y@meta.data
  # colnames(meta) <- c("barcode_id","orig.ident","nCount_RNA","nFeature_RNA","sample","disease","pathology","batch","patient","patient2","plaque","sex","age","pmi","seurat_clusters_ref","sample_id","test","RNA_snn_res.0.8")
  # y@meta.data <- meta
  # head(y@meta.data)
  # filter only the singlets
  # datasc_filter <- subset(y, subset = DF == "Singlet")
  saveRDS(object = y,file = paste0("../../out/object/WM_Nature2021_",x,"_SoupX.rds"))
})
