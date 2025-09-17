# martina wants to integrate the two tissue in one single object, in order to be able to score senescence

# LIBRARIES ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the WM data from Martina nature paper ---------------------------
# these are the samples reprocessed from the original analysis of martina (Nature 2021)
folder <- "../../out/object/"
file <- dir(folder) %>% str_subset(pattern = "WM_Nature2021_") %>%
  str_remove_all(pattern = ".rds") %>% 
  str_remove_all(pattern = "WM_Nature2021_|_SoupX") %>% 
  str_subset(pattern = "s25|s26|s27|s31",negate = T)

# read in the data
list_WM_new <- lapply(file,function(x){
  readRDS(paste0(folder,"WM_Nature2021_",x,"_SoupX.rds"))
}) %>% 
  setNames(file)

# test <- list_WM_new$s1
# name <- "s1"
# save the individual sample after uniforming the metadata
pmap(list(names(list_WM_new),list_WM_new),function(name,test){
  print(name)
  # make suer RNA is default
  DefaultAssay(test) <- "RNA"
  # extract the meta and add the sample name
  meta <- test@meta.data
  
  # uniform the meta
  meta_uniform <- meta[,c("sample_id","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters_ref")]
  
  # swap the meta
  test@meta.data <- meta_uniform
  
  # save the object
  obj_out <- paste0("../../data/individual_samples/run_all_SoupX/",name,".rds")
  saveRDS(test,file = obj_out)
})

# read in the WM new data -------------------------------------------------
# these are WM but they are not in the original analysis of martina (Nature 2021)
folder <- "../../out/object/"
file <- dir(folder) %>% str_subset(pattern = "WM_Nature2021_") %>%
  str_remove_all(pattern = ".rds") %>% 
  str_remove_all(pattern = "WM_Nature2021_|_SoupX") %>%
  str_subset(pattern = "s25|s26|s27|s31")

# read in the data
list_WM_new <- lapply(file,function(x){
  readRDS(paste0(folder,"WM_Nature2021_",x,"_SoupX.rds"))
}) %>% 
  setNames(file)

# test <- list_WM_new$s25
# name <- "s25"
# save the individual sample after uniforming the metadata
pmap(list(names(list_WM_new),list_WM_new),function(name,test){
  print(name)
  # make suer RNA is default
  DefaultAssay(test) <- "RNA"
  # extract the meta and add the sample name
  meta <- test@meta.data
  
  # uniform the meta
  meta_uniform <- meta[,c("sample_id","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters_ref")]
  
  # swap the meta
  test@meta.data <- meta_uniform
  
  # save the object
  obj_out <- paste0("../../data/individual_samples/run_all_SoupX/",name,".rds")
  saveRDS(test,file = obj_out)
})

# # read in the individual WM objects from Martina --------------------------
# folder <- "../../data/WM_data/WM_data_Abisnta_new/"
# file <- dir(folder) %>% str_subset(pattern = ".rds") %>% str_remove_all(pattern = ".rds")
# 
# # read in the data
# list_WM_new <- lapply(file,function(x){
#   readRDS(paste0(folder,x,".rds"))
# }) %>% 
#   setNames(file)
# 
# # harmonyze the meta and save the objects
# pmap(list(names(list_WM_new),list_WM_new),function(name,test){
#   print(name)
#   # make suer RNA is default
#   DefaultAssay(test) <- "RNA"
#   # extract the meta and add the sample name
#   meta <- test@meta.data %>% 
#     mutate(sample_id=name) %>% 
#     dplyr::rename(seurat_clusters_ref=seurat_clusters)
#   
#   # uniform the meta
#   meta_uniform <- meta[,c("sample_id","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters_ref")]
#   
#   # swap the meta
#   test@meta.data <- meta_uniform
#   
#   # save the object
#   obj_out <- paste0("../../data/individual_samples/",name,".rds")
#   saveRDS(test,file = obj_out)
# })

# read in the data for the cortex datasets --------------------------------
CX_ansinta <- readRDS("../../out/object/ManualClean/data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony.rds")
# DimPlot(CX_ansinta,raster = T,group.by = "seurat_clusters",label = T)
# split the object into the individula samples
# meta_cx <- CX_ansinta@meta.data %>% 
#   rownames_to_column("barcode")
# 
# table(meta_cx$official_id,meta_cx$orig.ident.cca)

# split the object into individula samples
list_CX_absinta <- SplitObject(CX_ansinta,split.by = "official_id")

# save the individual sample after uniforming the metadata
pmap(list(names(list_CX_absinta),list_CX_absinta),function(name,test){
  print(name)
  # make suer RNA is default
  DefaultAssay(test) <- "RNA"
  # extract the meta and add the sample name
  meta <- test@meta.data %>% 
    mutate(sample_id=name) %>% 
    dplyr::rename(seurat_clusters_ref=seurat_clusters,
                  percent.mt = percent.mt.harmony)
  
  # uniform the meta
  meta_uniform <- meta[,c("sample_id","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters_ref")]
  
  # swap the meta
  test@meta.data <- meta_uniform
  
  # save the object
  obj_out <- paste0("../../data/individual_samples/run_all_SoupX/",name,".rds")
  saveRDS(test,file = obj_out)
})

# chek the total dataset --------------------------------------------------
df_sample <- data.frame(official_id = dir("../../data/individual_samples/run_all_SoupX/") %>% 
                          str_remove_all(".rds")) %>% 
  mutate(present = 1)

df_tot_sample <- read_csv("../../data/LUT_sample_WM_CX.csv")

df_meta_sample1 <- left_join(df_tot_sample,df_sample)

df_meta_sample1 %>% 
  group_by(present,pathology_class) %>% 
  summarise(n=n())

df_meta_sample1 %>% 
  filter(is.na(present))

df_meta_sample1 %>%
  filter(origin == "wm") %>% 
  group_by(pathology_class,official_id) %>% 
  summarise(n=n()) %>% 
  print(n=30)
