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
data.combined <- readRDS("../../out/object/revision/100_subcluster_OLIGO.rds")
Idents(data.combined) <- "seurat_clusters"

# read in the original dataset
WM_ansinta_nature <- readRDS("../../data/WM_data/WM_data_Absinta_Nature2021/all20_integrated_clean_metadata.rds")

# plots -------------------------------------------------------------------
# general UMAP with new clustering
DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.2", split.by = "origin",label = T,raster = T,ncol=5)
DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.2", split.by = "origin_alt",label = T,raster = T,ncol=5)
# DimPlot(data.combined, reduction = "umap", split.by = "origin",label = T,raster = T,ncol=5)
# general UMAP with former clustering
DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters_refCXWM",split.by = "origin",label = T,raster = T,ncol=5)

# general UMAP with former clustering
DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters_ref",split.by = "origin",label = T,raster = T,ncol=5)

# subset only the original data from Martina
test <- subset(data.combined,subset = origin_alt == "wm")
test$seurat_clusters_ref <- test$seurat_clusters_ref %>%
  fct_relevel(as.character(0:17))

# plot in comparison with the new embedding
(DimPlot(WM_ansinta_nature,group.by = "seurat_clusters",raster=T,label = T)+ggtitle("original UMAP"))+
  (DimPlot(test,group.by = "seurat_clusters_ref",split.by = "origin_alt",raster=T,label = T) + ggtitle("new UMAP"))

DimPlot(test,group.by = "seurat_clusters_ref",split.by = "seurat_clusters_ref",raster=F,label = T,ncol=4) + ggtitle("new UMAP")

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

wrap_plots(list_plot)
# ggsave("../../out/plot/manualClean/UMAPCluster_resolutions.pdf",width = 25,height = 15)

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

# make a shortlist of the markers
shortlist_features_list_short <- list(
  IMMUNE = c("CX3CR1","P2RY12","C3","CSF1R", "CD74","C1QB"),
  B_CELLS = c("IGHG1", "CD38"),
  T_CELLS =  c("SKAP1", "CD8A", "CD2"),
  OLIGO = c("MOG","MBP","MAG"),
  OPC = c("NLGN4X","OLIG1","OLIG2"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1"),
  NEURONS = c("CUX2","SYP", "NEFL","SYT1"),
  VAS = c("VWF","FLT1","CLDN5","PDGFRB"),
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
                       group.by = "RNA_snn_res.0.2") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.2")+
  theme(strip.text = element_text(angle = 90))
# ggsave(plot=test_long01,"../../out/plot/manualClean/DotplotLong_harmonySkipIntegration_AllSoupX_00500_07000_05_res0.1.pdf",width = 30,height = 6)

# score the signatures, both long and short
data.combined <- Seurat::AddModuleScore(data.combined,
                                        features = shortlist_features_list_short,
                                        name = "_score_short")

data.combined <- Seurat::AddModuleScore(data.combined,
                                        features = shortlist_features_list_long,
                                        name = "_score_long")

df_rename_short <- data.frame(names = data.combined@meta.data %>%
                                colnames() %>%
                                str_subset("_score_short"),
                              rename = paste0("scoreShort_",names(shortlist_features_list_short)))

df_rename_long <- data.frame(names = data.combined@meta.data %>%
                               colnames() %>%
                               str_subset("_score_long"),
                             rename = paste0("scoreLong_",names(shortlist_features_list_long)))

lookup_short <- df_rename_short$names
names(lookup_short) <- df_rename_short$rename

lookup_long <- df_rename_long$names
names(lookup_long) <- df_rename_long$rename

# rename the columns
data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup_short))
data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup_long))

# plot the scores from AddModuleScore
list_plot_02_short <- lapply(df_rename_short$rename,function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_short)
# ggsave("../../out/plot/manualClean/UMAPCluster_ExpertAnnotaiton_short.pdf",width = 22,height = 15)

list_plot_02_long <- lapply(df_rename_long$rename,function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_long)
# ggsave("../../out/plot/manualClean/UMAPCluster_ExpertAnnotaiton_long.pdf",width = 22,height = 15)

# same as above but as violin plot
list_plot <- lapply(df_rename_long$rename, function(x){ 
  test <- VlnPlot(object = data.combined,features = x, group.by = "RNA_snn_res.0.2",raster = T)
  return(test)
})

# make it a dataframe
x <- list_plot[[1]]
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

# how many cells per ident
df_violin %>%
  group_by(ident,feature) %>%
  summarise(n = n()) %>%
  filter(feature %in% c("score_ASTRO"))

# plot at maximum 1000 cells per group

df_violin %>% 
  group_by(ident,feature) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = ident,values_from = n) %>%
  print(n=11)

# fi the sample size of the cluster is below 1000 cells use them all, else randomly select a subset 0f 1000 cells
test <- df_violin %>%
  group_by(ident,feature) %>%
  nest() %>%
  mutate(n_cell = map(data,function(x){
    dim(x)[1] 
  }) %>% unlist())

set.seed(123)
df_plot_violin <- bind_rows(test %>%
            filter(n_cell>1000) %>%
            mutate(data_subset = map(data,function(x){
              x %>%
                sample_n(size = 400,replace = F)
            })),
          
          test %>%
            filter(n_cell<1000) %>%
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
# ggsave("../../out/plot/manualClean/ViolinCluster_ExpertAnnotaiton_res0.2.pdf",width = 7,height = 21)


# marker clusters ---------------------------------------------------------
# set the Idents
Idents(data.combined) <- "RNA_snn_res.0.2"

res_test_res0.2 <- SeuratWrappers::RunPrestoAll(data.combined,only.pos = T,logfc.threshold = 0.25,min.pct = 0.1)
# save the table
write_tsv(res_test_res0.2,"../../out/table/revision/100_subcluster_OLIGO_res_test_res0.2.tsv")

res_test_res0.2 %>%
  filter(cluster == 2)

res_test_res0.2 %>%
  filter(cluster == 0)

# run enrichr with the list of genes in the module
# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
# filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "GO"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023")

# query -------------------------------------------------------------------
list_genes_UP <- list(clust_2_UP = res_test_res0.2 %>%
                        filter(cluster == 2) %>% pull(gene))

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_enrichr_UP <- lapply(list_genes_UP,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>%
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_enrichr_UP %>%
  write_tsv("../../out/table/revision/100_subcluster_OLIGO_EnrichR_clust2_res0.2.tsv")


# list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")
plot_list_UP <- list_enrichr_UP %>%
  split(f = .$comparison)

# library(scales)
list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
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

wrap_plots(list_plot_UP,nrow = 1)
# ggsave("../../out/image/enrichR_DE_pseudobulk_ASTRO_ctrl_refCX_UP_shr_filterExp.pdf",width = 6,height = 12,limitsize = FALSE)