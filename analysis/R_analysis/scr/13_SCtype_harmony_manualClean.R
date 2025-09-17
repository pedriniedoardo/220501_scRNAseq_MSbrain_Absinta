# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggraph)
library(igraph)
library(data.tree)
library(HGNChelper)
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R",destfile = "scr/R/gene_sets_prepare.R")
source("scr/gene_sets_prepare.R")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R",destfile = "scr/R/sctype_score_.R")
source("scr/sctype_score_.R")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R",destfile = "scr/R/auto_detect_tissue_type.R")
source("scr/auto_detect_tissue_type.R")

# download the database file
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",destfile = "data/ScTypeDB_full.xlsx")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",destfile = "data/ScTypeDB_short.xlsx")

# read in the final object ------------------------------------------------
scobj <- readRDS("../../out/object/ManualClean/data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony.rds")
# sobj_test_cluster <- FindClusters(scobj,resolution = 0.3)

# data.combined <- readRDS("../out_large/scRNAseq_analysis/object/data.combined_fix_filter_norm_DoubletSinglet_integrated_5K.rds")

# confirm the identity of the object
DimPlot(scobj,label = T,raster = T)

# DimPlot(sobj_test_cluster,label = T)

# confirm the object has a scaled matrix in it. this is the one that will bu used for the scoring of the cell types
# data.combined@assays$RNA@scale.data[1:10,1:10]
# dim(data.combined@assays$RNA@scale.data)
scobj@assays$RNA@scale.data[1:10,1:10]
dim(scobj@assays$RNA@scale.data)

# read in the database file
db_ <- "../../data/ScTypeDB_full.xlsx"
tissue <- "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)
str(gs_list)

# Finally, let's assign cell types to each cluster:
# get cell-type by cell matrix
es.max <- sctype_score(scRNAseqData = scobj[["RNA"]]@scale.data, scaled = TRUE, 
                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

es.max[1:10,1:10]

# wrangling ---------------------------------------------------------------
# after running the annotation in theory all the cells can have a specific annotation
df_score <- es.max %>% 
  data.frame() %>% 
  rownames_to_column("annotation") %>% 
  pivot_longer(names_to = "barcode",values_to = "score",-annotation) %>% 
  mutate(barcode = str_replace(barcode,pattern = "\\.",replacement = "-"))

df_score
# add the new annotation to the original object and show the umap
df_meta <- scobj@meta.data %>% 
  rownames_to_column("barcode")

head(df_meta)

# for each cluster in the dataset sum the scores per cell for each cluster
df_score_cluster <- lapply(unique(df_meta$seurat_clusters),function(cl){
  # pull the cells belonging to the clster
  id <- df_meta %>% 
    filter(seurat_clusters %in% cl) %>% 
    pull(barcode)
  
  # sum the score from eah cell for each category
  df_tot <- df_score %>% 
    filter(barcode %in% id) %>% 
    # add the cluster information
    mutate(cluster = cl) %>% 
    # add the total number of cells
    # mutate(n_cell = length(id)) %>% 
    group_by(annotation,cluster) %>% 
    summarise(sum_score = sum(score),
              n_cell = n())
  
  return(df_tot)
  
}) %>% 
  bind_rows()

df_score_cluster

# pull the top score per cluster
df_score_cluster_top <- df_score_cluster %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = sum_score) %>% 
  # set low-confident (low ScType score) clusters to "unknown"
  mutate(annotation_confident = case_when(sum_score < n_cell/4 ~"Unknown",
                                          T~annotation)) %>% 
  ungroup()

df_score_cluster_top

# add the annotation to the original dataset
df_meta_full <- left_join(df_meta,df_score_cluster_top,c("seurat_clusters"="cluster"))
head(df_meta_full)

# add it back to the object
scobj@meta.data <- df_meta_full %>% 
  column_to_rownames("barcode")

head(scobj@meta.data)

# plotting ----------------------------------------------------------------
DimPlot(scobj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'annotation',raster = T)
ggsave("../../out/image/ManualClean/UMAP_data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_SCType.pdf",width = 10,height = 7)

# save the annotated object -----------------------------------------------
saveRDS(scobj,"../../out/object/ManualClean/data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_SCType.rds")

