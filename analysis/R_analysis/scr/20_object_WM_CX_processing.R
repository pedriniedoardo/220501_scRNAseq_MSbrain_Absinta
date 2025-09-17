# martina wants to integrate the two tissue in one single object, in order to be able to score senescence

# LIBRARIES ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the WM data from Martina nature paper ---------------------------
WM_ansinta_nature <- readRDS("../../data/WM_data/WM_data_Absinta_Nature2021/all20_integrated_clean_metadata.rds")

# read in the corrected LUT for the sample names
LUT_sampleNature2021<- read_csv("../../data/LUT_WM_Nature2021.csv") %>% 
  mutate(sample_seurat = as.character(sample_seurat))

# split the object into the individula samples
meta_nature <- WM_ansinta_nature@meta.data %>% 
  rownames_to_column("barcode") %>% 
  left_join(LUT_sampleNature2021,by = c("sample"="sample_seurat"))

table(meta_nature$sample_id,meta_nature$pathology)

# add the new meta to the original object
WM_ansinta_nature$sample_id <- meta_nature$sample_id

# split the object into individula samples
list_WM_nature <- SplitObject(WM_ansinta_nature,split.by = "sample_id")
# test<-list_WM_nature$s1
# save the individual sample after uniforming the metadata
pmap(list(names(list_WM_nature),list_WM_nature),function(name,test){
  print(name)
  # make suer RNA is default
  DefaultAssay(test) <- "RNA"
  # extract the meta and add the sample name
  meta <- test@meta.data %>% 
    mutate(sample_id=name) %>% 
    dplyr::rename(seurat_clusters_ref=seurat_clusters)
  
  # uniform the meta
  meta_uniform <- meta[,c("sample_id","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters_ref")]
  
  # swap the meta
  test@meta.data <- meta_uniform
  
  # save the object
  obj_out <- paste0("../../data/individual_samples/",name,".rds")
  saveRDS(test,file = obj_out)
})

# read in the individual WM objects from Martina --------------------------
folder <- "../../data/WM_data/WM_data_Abisnta_new/"
file <- dir(folder) %>% str_subset(pattern = ".rds") %>% str_remove_all(pattern = ".rds")

# read in the data
list_WM_new <- lapply(file,function(x){
  readRDS(paste0(folder,x,".rds"))
}) %>% 
  setNames(file)

# harmonyze the meta and save the objects
pmap(list(names(list_WM_new),list_WM_new),function(name,test){
  print(name)
  # make suer RNA is default
  DefaultAssay(test) <- "RNA"
  # extract the meta and add the sample name
  meta <- test@meta.data %>% 
    mutate(sample_id=name) %>% 
    dplyr::rename(seurat_clusters_ref=seurat_clusters)
  
  # uniform the meta
  meta_uniform <- meta[,c("sample_id","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters_ref")]
  
  # swap the meta
  test@meta.data <- meta_uniform
  
  # save the object
  obj_out <- paste0("../../data/individual_samples/",name,".rds")
  saveRDS(test,file = obj_out)
})

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
  obj_out <- paste0("../../data/individual_samples/",name,".rds")
  saveRDS(test,file = obj_out)
})

# chek the total dataset --------------------------------------------------
df_sample <- data.frame(official_id = dir("../../data/individual_samples/") %>% 
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
  