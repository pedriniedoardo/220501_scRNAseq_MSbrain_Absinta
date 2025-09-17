# AIM ---------------------------------------------------------------------
# try the subclustering of the OLIGO cells

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
library(scales)
library(cowplot)
library(SeuratWrappers)
library(enrichR)

# read the data -----------------------------------------------------------
# read in the dataset
data.combined <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegAllSoupX_expertAnno.rds")

# save the plot
DimPlot(data.combined, group.by = "expertAnno.l1",raster=T,label = T)
dim(data.combined)

# subset the cells
data.combined_subset <- subset(data.combined,subset = expertAnno.l1 == "OLIG")

# read in the data --------------------------------------------------------

# save the sperse matrix
total_sm <- data.combined_subset@assays$RNA@counts

# check the total_sm
dim(total_sm)
# confirmt he total number of cells

# check a sample from the dataset
total_sm[1:10,50:70]

# I need to create a single object to add the cell cycle scoring and other metadata
# I want to keep all the cells

# I need to create a single object to add the cell cycle scoring and other metadata
sobj_total <- CreateSeuratObject(counts = total_sm,
                                 meta.data = data.combined_subset@meta.data %>%
                                   dplyr::rename(seurat_clusters_old = seurat_clusters),
                                 project = "CX_WM_OLIGO",
                                 # keep all the filtered cells
                                 min.cells = 0, min.features = 0)

# add the cell cycle analysis
# DefaultAssay(sobj_total) <- "RNA"
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes, g2m.features = g2m.genes)
# sobj_total$percent.mt <- PercentageFeatureSet(sobj_total, pattern = "^MT-")
# sobj_total$percent.ribo <- PercentageFeatureSet(sobj_total, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

# # add all the original metadata 
# meta_new <- sobj_total@meta.data %>% 
#   rownames_to_column("barcode")
# 
# meta_new_total <- left_join(meta_new,df_meta,"barcode",suffix=c(".harmony",".cca")) %>% 
#   left_join(LUT,by = c("orig.ident"="official_id")) %>% 
#   column_to_rownames("barcode")
# 
# # update the meta
# sobj_total@meta.data <- meta_new_total

# pull all the variable to scale
sobj_total@assays$RNA@scale.data
# all.genes <- rownames(sobj_total)

# rescale the data for regressing out the sources of variation do not scale all the genes. if needed scale them before the heatmap call
sobj_total <- sobj_total %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  # I can scale the missing features afterwards
  ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score","G2M.Score","origin","facility"), verbose = T) %>% 
  # ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T,features = all.genes) %>% 
  RunPCA(npcs = 30, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2) %>%
  identity()

# check the object before harmony integration
# save the plot
DimPlot(sobj_total, group.by = "orig.ident",raster=T,label = T)
DimPlot(sobj_total, group.by = "pathology_class",raster=T,label = T)
DimPlot(sobj_total, group.by = "origin",raster=T,label = T)

saveRDS(sobj_total,"../../out/object/ManualClean/data.combined_WM_CX_OLIGO_SkipIntegration.rds")

table(sobj_total@meta.data$pathology_class)

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

DimPlot(sobj_total_h, group.by = "orig.ident",raster=T,label = T)
DimPlot(sobj_total_h, group.by = "pathology_class",raster=T,label = T)
DimPlot(sobj_total_h, group.by = "origin",raster=T,label = T)
DimPlot(sobj_total_h,split.by = "origin", group.by = "seurat_clusters",raster=T,label = T)

# saveRDS(sobj_total_h,"../out_large/scRNAseq_analysis/object/sobj_total_h_fix_filter_norm_doublet_harmony_5K.rds")
saveRDS(sobj_total_h,"../../out/object/ManualClean/data.combined_WM_CX_OLIGO_harmonySkipIntegration.rds")

# run DE per cluster ------------------------------------------------------
sobj_total_h.markers <- RunPrestoAll(sobj_total_h, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_tsv(sobj_total_h.markers,"../../out/table/ManualClean/FindAllMarkers_data.combined_WM_CX_OLIGO_harmonySkipInteg_PosOnly.tsv")

# run enrichR per cluster -------------------------------------------------
dbs <- listEnrichrDbs()

#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Tabula"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Azim"))

dbs_db <- c("KEGG_2021_Human",
            "MSigDB_Hallmark_2020",
            "Reactome_2016",
            "HDSigDB_Human_2021","Tabula_Sapiens","Azimuth_Cell_Types_2021")

# GENE SELECTION ----------------------------------------------------------
#
list_res_tot_UP <- sobj_total_h.markers %>%
  split(f = .$cluster) %>%
  lapply(function(x){
    x %>%
      filter(pct.1 - pct.2 > 0.15,
             avg_log2FC > 1) %>%
      pull(gene)
  })

# list_res_tot_UP <- sobj_total_h.markers %>%
#   split(f = .$cluster) %>%
#   lapply(function(x){
#     x %>%
#       filter(avg_log2FC > 1) %>%
#       pull(gene)
#   })

# query -------------------------------------------------------------------
# seelct only the clusters with more than 10 genes as degs
list_res_tot_UP_filter <- list_res_tot_UP[lengths(list_res_tot_UP)>5]

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_UP <- lapply(list_res_tot_UP_filter,function(x){
  genes <- x
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_UP %>%
  write_tsv("../../out/table/ManualClean/enrichR_FindAllMarkers_WM_CX_OLIGO_harmonySkipInteg_PosOnly.tsv")

# explore a specific term
list_UP %>%
  filter(annotation %in% c("HDSigDB_Human_2021"),
         comparison %in% c(1))
# filter(str_detect(Term,pattern = "OLIGO")) %>% 
#   dplyr::select(annotation,Term,Genes)

# list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")
plot_list_UP <- list_UP %>%
  split(f = .$comparison)

# library(scales)
list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_UP)
ggsave("../../out/image/ManualClean/enrichR_FindAllMarkers_WM_CX_OLIGO_harmonySkipInteg_PosOnly.pdf",width = 20,height = 35,limitsize = FALSE)
