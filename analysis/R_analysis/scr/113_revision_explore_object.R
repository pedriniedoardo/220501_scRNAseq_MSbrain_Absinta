# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration.
# this integration was run by skipping the seurat integration (by merging the matrices) and running Harmony
# to be fixed

# libraries ---------------------------------------------------------------
library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(SeuratWrappers)
library(dittoSeq)
library(clustree)
library(pals)
library(patchwork)
library(ComplexHeatmap)
library(enrichR)


# read in the data --------------------------------------------------------
# read in the original dataset
WM_ansinta_nature <- readRDS("../../data/WM_data/WM_data_Absinta_Nature2021/all20_integrated_clean_metadata.rds")

# loop processing ---------------------------------------------------------
# define the resolution parameter
res_id <- c(1000,2000,3000,4000)
# res <- 8000

# loop the processing
lapply(res_id, function(res){
  # read in the data --------------------------------------------------------
  data.combined <- readRDS(paste0("../../out/object/revision/112_WMCX_ManualClean2_harmonySkipIntegration_AllSoupX_",res,".rds"))
  Idents(data.combined) <- "seurat_clusters"
  
  # plots -------------------------------------------------------------------
  # general UMAP with new clustering
  # DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.2", split.by = "origin",label = T,raster = T,ncol=5)
  # DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.2", split.by = "origin_alt",label = T,raster = T,ncol=5)
  # DimPlot(data.combined, reduction = "umap", split.by = "origin",label = T,raster = T,ncol=5)
  
  # plot the data with the previous defiend version of the dataset
  # DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters.old",split.by = "seurat_clusters.old",label = T,raster = T,ncol=4)
  DimPlot(data.combined, reduction = "umap", group.by = "expertAnno.l2.old",split.by = "origin",label = T,raster = T,ncol=4) +plot_annotation(paste0("HGV test ",res))
  ggsave(paste0("../../out/image/revision/113_",res,"_expertAnno.l2.old.pdf"),width = 12,height = 6)
  
  DimPlot(data.combined, reduction = "umap", group.by = "RNA_snn_res.0.3",split.by = "origin",label = T,raster = T,ncol=4) +plot_annotation(paste0("HGV test ",res))
  ggsave(paste0("../../out/image/revision/113_",res,"_RNA_snn_res.0.3.pdf"),width = 12,height = 6)
  
  # general UMAP with former clustering
  # DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters_ref",split.by = "origin",label = T,raster = T,ncol=5)
  
  # subset only the original data from Martina
  test <- subset(data.combined,subset = origin_alt == "wm")
  test$seurat_clusters_ref <- test$seurat_clusters_ref %>%
    fct_relevel(as.character(0:17))
  
  # plot in comparison with the new embedding
  (DimPlot(WM_ansinta_nature,group.by = "seurat_clusters",raster=T,label = T)+ggtitle("original UMAP"))+
    (DimPlot(test,group.by = "seurat_clusters_ref",split.by = "origin_alt",raster=T,label = T) + ggtitle("new UMAP"))+plot_annotation(paste0("HGV test ",res))
  ggsave(paste0("../../out/image/revision/113_",res,"_cluster_Nature2021.pdf"),width = 12,height = 6)
  
  DimPlot(test,group.by = "seurat_clusters_ref",split.by = "seurat_clusters_ref",raster=T,label = T,ncol=4) + ggtitle(paste0("New HGV test ",res))
  ggsave(paste0("../../out/image/revision/113_",res,"_cluster_Nature2021_split.pdf"),width = 15,height = 15)
  
  # plot the UMAP with all the resolutions runs
  id_resolution <- str_subset(colnames(data.combined@meta.data),pattern = "RNA_snn_res") %>%
    sort()
  
  list_plot <- lapply(id_resolution,function(x){
    plot <- DimPlot(data.combined,
                    reduction = "umap",
                    group.by = x,
                    label = T,
                    raster = T)
    return(plot)
  })
  
  wrap_plots(list_plot) + plot_annotation(paste0("HGV test ",res))
  ggsave(paste0("../../out/image/revision/113_",res,"_cluster_res_split.pdf"),width = 25,height = 15)
  
  # martina also suggested to reomve the NPC
  shortlist_features_list_long <- list(
    IMMUNE = c("LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
    B_CELLS = c("IGHG1", "CD38"),
    T_CELLS =  c("SKAP1", "CD8A", "CD2"),
    OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG"),
    ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
    NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
    ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
    PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"),
    SCHWANN = c("PMP22","MPZ","PRX"),
    EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
    STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
  )
  
  # plot the shortlisted feature per cluster
  # notice that this is done only on the subset of the young (control) cells
  test_long01 <- DotPlot(data.combined,
                         features = shortlist_features_list_long,
                         dot.scale = 8,
                         cluster.idents = T,
                         group.by = "RNA_snn_res.0.3") +
    RotatedAxis() +
    labs(title = paste0("HGV test ",res," RNA_snn_res.0.3"))+
    theme(strip.text = element_text(angle = 90))
  ggsave(plot = test_long01, paste0("../../out/image/revision/113_",res,"_dotplot_res.0.3.pdf"),width = 25,height = 7)
  
  # score the signatures, both long and short
  data.combined <- Seurat::AddModuleScore(data.combined,
                                          features = shortlist_features_list_long,
                                          name = "_score_long")
  
  df_rename_long <- data.frame(names = data.combined@meta.data %>%
                                 colnames() %>%
                                 str_subset("_score_long"),
                               rename = paste0("scoreLong_",names(shortlist_features_list_long)))
  
  lookup_long <- df_rename_long$names
  names(lookup_long) <- df_rename_long$rename
  
  # rename the columns
  data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup_long))
  
  list_plot_02_long <- lapply(df_rename_long$rename,function(x){
    plot <- FeaturePlot(data.combined,features = x,order = T,
                        reduction = "umap",
                        raster = T) + scale_color_viridis_c(option = "turbo")
    return(plot)
  })
  
  wrap_plots(list_plot_02_long) + plot_annotation(paste0("HGV test ",res))
  ggsave(paste0("../../out/image/revision/113_",res,"_signatureScore.pdf"),width = 22,height = 15)
  
  # same as above but as violin plot
  list_plot <- lapply(df_rename_long$rename, function(x){ 
    test <- VlnPlot(object = data.combined,features = x, group.by = "RNA_snn_res.0.3",raster = T)
    return(test)
  })
  
  # make it a dataframe
  # x <- list_plot[[1]]
  df_violin <- lapply(list_plot,function(x){ 
    df <- x[[1]]$data 
    
    # extract the name of the gene 
    feature <- colnames(df)[1] 
    
    df %>% 
      mutate(feature = feature) %>% 
      setNames(c("value","ident","feature")) 
  }) %>% 
    bind_rows()
  
  head(df_violin) 
  
  # fi the sample size of the cluster is below 1000 cells use them all, else randomly select a subset 0f 1000 cells
  test_violin <- df_violin %>%
    group_by(ident,feature) %>%
    nest() %>%
    mutate(n_cell = map(data,function(x){
      dim(x)[1] 
    }) %>% unlist())
  
  set.seed(123)
  df_plot_violin <- bind_rows(test_violin %>%
                                dplyr::filter(n_cell>500) %>%
                                mutate(data_subset = map(data,function(x){
                                  x %>%
                                    sample_n(size = 500,replace = F)
                                })),
                              
                              test_violin %>%
                                dplyr::filter(n_cell<500) %>%
                                mutate(data_subset = data)) %>%
    unnest(data_subset) %>%
    select(ident,feature,value) %>%
    ungroup()
  
  df_plot_violin_summary <- df_plot_violin %>%
    group_by(feature) %>%
    summarise(med_score = median(value))
  
  df_plot_violin %>%
    ggplot(aes(y = ident, x = value)) + 
    geom_violin(scale = "width")+ 
    #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) + 
    geom_point(position = position_jitter(width = 0.2),alpha = 0.05,size = 0.5) + 
    facet_wrap(~feature,nrow = 1,scales = "free") + 
    theme_bw() + 
    geom_vline(data = df_plot_violin_summary,aes(xintercept = med_score),linetype="dashed",col="red") +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(hjust = 1,angle = 45))
  ggsave(paste0("../../out/image/revision/113_",res,"_signatureScore_violin.pdf"),width = 22,height = 15)
  
  
  # marker clusters ---------------------------------------------------------
  # set the Idents
  Idents(data.combined) <- "RNA_snn_res.0.3"
  
  res_test_res0.3 <- SeuratWrappers::RunPrestoAll(data.combined,only.pos = T,logfc.threshold = 0.25,min.pct = 0.1)
  # save the table
  write_tsv(res_test_res0.3,paste0("../../out/table/revision/113_",res,"_AllMarkers_res.0.3.tsv"))
})
