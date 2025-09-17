# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration.
# this integration was run by skipping the seurat integration (by merging the matrices) and running Harmony

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

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/revision/118_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000.rds")
Idents(data.combined) <- "seurat_clusters"

# data.combined@meta.data %>%
#   filter(RNA_snn_res.0.1 %in% c(8,9)) %>%
#   group_by(RNA_snn_res.0.1,cluster_whitelist) %>%
#   summarise(n = n()) %>%
#   arrange(RNA_snn_res.0.1,desc(n))

DimPlot(data.combined,raster=T,label=T,split.by ="origin",group.by = "expertAnno.l1.old")
DimPlot(data.combined,raster=T,label=T,split.by ="origin",group.by = "expertAnno.l2.old")
DimPlot(data.combined,raster=T,label=T,split.by ="origin",group.by = "RNA_snn_res.0.2")
DimPlot(data.combined,raster=T,label=T,split.by ="origin",group.by = "seurat_clusters.old")
DimPlot(data.combined,raster=T,label=T,split.by ="origin_alt",group.by = "RNA_snn_res.0.2")
DimPlot(data.combined,raster=T,label=T,split.by ="origin_alt",group.by = "seurat_clusters_ref")

# plots -------------------------------------------------------------------
# plot the tree of the cluster dependencies. this will justify the choice of the resolution, not too granula for the moment.
# library(clustree)
clustree::clustree(data.combined@meta.data[,grep("RNA_snn_res", colnames(data.combined@meta.data))],
                   prefix = "RNA_snn_res.")
ggsave("../../out/image/revision/119_UMAPCluster_tree.pdf",width = 10,height = 10)

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
ggsave("../../out/image/revision/119_UMAPCluster_resolutions.pdf",width = 25,height = 15)

# save the plot with the specific resolutiun selected. use the res 0.3
DimPlot(data.combined,
        reduction = "umap",
        group.by = "RNA_snn_res.0.3",
        label = T,
        raster = T)
ggsave("../../out/image/revision/119_UMAPCluster_resolutions_0.3.pdf",width = 6,height = 5)

DimPlot(data.combined,
        reduction = "umap",
        group.by = "RNA_snn_res.0.3",split.by = "origin",
        label = T,
        raster = T)
ggsave("../../out/image/revision/119_UMAPCluster_resolutions_0.3_origin.pdf",width = 11,height = 5)

DimPlot(data.combined,
        reduction = "umap",
        group.by = "RNA_snn_res.0.1",
        label = T,
        raster = T)
ggsave("../../out/image/revision/119_UMAPCluster_resolutions_0.1.pdf",width = 6,height = 5)

DimPlot(data.combined,
        reduction = "umap",
        group.by = "RNA_snn_res.0.1",split.by = "origin",
        label = T,
        raster = T)
ggsave("../../out/image/revision/119_UMAPCluster_resolutions_0.1_origin.pdf",width = 11,height = 5)

# save the panel martina uses for the annotation
shortlist_features_list_long <- list(
  IMMUNE = c("LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  B_CELLS = c("IGHG1", "CD38"),
  T_CELLS =  c("SKAP1", "CD8A", "CD2"),
  OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
  NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
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

# make a tailored one following martina's reccomendation
shortlist_features_list_short2 <- list(
  IMMUNE = c("CX3CR1","P2RY12","CSF1R","RUNX1","SKAP1","PTPRC"),
  OLIGO = c("PLP1","MBP","MOG"),
  OPC = c("SOX6","PDGFRA","NLGN4X"),
  ASTRO = c("AQP4","SLC1A2","GFAP"),
  VAS = c("CLDN5","FLT1"),
  NEURONS = c("SYT1","GAD2","RORB","SATB2","CUX2","NRGN","SLC17A7"),
  EPENDYMA = c("CFAP299","DNAH7","DNAH9")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long01 <- DotPlot(data.combined,
                       features = shortlist_features_list_long,
                       dot.scale = 8,
                       cluster.idents = T,
                       group.by = "RNA_snn_res.0.3") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.3")+
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_long01,"../../out/image/revision/119_DotplotLong_res0.3.pdf",width = 30,height = 6)

# plot the shortlist
test_short04 <- DotPlot(data.combined,
                        features = shortlist_features_list_short,
                        dot.scale = 8,
                        cluster.idents = T,
                        group.by = "RNA_snn_res.0.3") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.3") +
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_short04,"../../out/image/revision/119_DotplotShort_res0.3.pdf",width = 15,height = 6)

# plot the shortlist 2
test_short05 <- DotPlot(data.combined,
                        features = shortlist_features_list_short2,
                        dot.scale = 8,
                        cluster.idents = T,
                        group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1") +
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_short05,"../../out/image/revision/119_DotplotShort2_res0.1.pdf",width = 10,height = 5)

test_short06 <- DotPlot(data.combined,
                        features = shortlist_features_list_short2,
                        dot.scale = 8,
                        cluster.idents = T,
                        group.by = "RNA_snn_res.0.3") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.3") +
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_short06,"../../out/image/revision/119_DotplotShort2_res0.3.pdf",width = 10,height = 6)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_short062 <- DotPlot(data.combined,
                        features = shortlist_features_list_short2,
                        dot.scale = 8,
                        cluster.idents = T,
                        group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1") +
  theme(strip.text = element_text(angle = 90))
# ggsave("../../out/image/ManualClean/Dotplot_data.combined_WM_CX_harmonySkipIntegAllSoupX_tailored.pdf",width = 12,height = 4)
# force the order suggested by matrina
df_test <- lapply(shortlist_features_list_short2,function(x){
  test_short062$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

glimpse(df_test)

df_test %>%
  # force the order
  mutate(id = factor(id,levels = c(3,11,0,6,2,12,13,9,4,5,7,8,1,10))) %>% 
  mutate(cell_type = factor(cell_type,levels = c("IMMUNE","OLIGO","OPC","ASTRO","VAS","NEURONS","EPENDYMA"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")
ggsave("../../out/image/revision/119_DotplotShort2_res0.1_tailoredReorder.pdf",width = 11,height = 4)

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
ggsave("../../out/image/revision/119_UMAPCluster_ExpertAnnotaiton_short.pdf",width = 22,height = 15)

list_plot_02_long <- lapply(df_rename_long$rename,function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_long)
ggsave("../../out/image/revision/119_UMAPCluster_ExpertAnnotaiton_long.pdf",width = 22,height = 15)

# same as above but as violin plot
list_plot <- lapply(df_rename_long$rename, function(x){ 
  test <- VlnPlot(object = data.combined,features = x, group.by = "RNA_snn_res.0.1",raster = T)
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

# plot at maximum 1000 cells per group
test_violin_score <- df_violin %>%
  group_by(ident,feature) %>%
  nest() %>%
  mutate(n_cell = map(data,function(x){
    dim(x)[1] 
  }) %>% unlist())

# plot at maximum 500 cells per group
set.seed(123)
df_plot_violin <- bind_rows(test_violin_score %>%
                              dplyr::filter(n_cell>500) %>%
                              mutate(data_subset = map(data,function(x){
                                x %>%
                                  sample_n(size = 500,replace = F)
                              })),
                            
                            test_violin_score %>%
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
ggsave("../../out/image/revision/119_ViolinCluster_ExpertAnnotaiton_res0.1.pdf",width = 21,height = 7)

# use a similar approach to score potential technical clusters
# plot the scores from AddModuleScore
list_plot_technical <- lapply(c("nFeature_RNA","percent.mt","percent.ribo","percent.globin"),function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_technical)
ggsave("../../out/image/revision/119_UMAPCluster_technical.pdf",width = 10,height = 8)

list_plot_technical <- lapply(c("nFeature_RNA","percent.mt","percent.ribo","percent.globin"), function(x){ 
  test <- VlnPlot(object = data.combined,features = x, group.by = "RNA_snn_res.0.1",raster = T)
  return(test)
})

# make it a dataframe
# x <- list_plot_technical[[1]]
df_violin_technical <- lapply(list_plot_technical,function(x){ 
  df <- x[[1]]$data 
  
  # extract the name of the gene 
  feature <- colnames(df)[1] 
  
  df %>% 
    mutate(feature = feature) %>% 
    setNames(c("value","ident","feature")) 
}) %>% 
  bind_rows()

head(df_violin_technical) 

test_violin <- df_violin_technical %>%
  group_by(ident,feature) %>%
  nest() %>%
  mutate(n_cell = map(data,function(x){
    dim(x)[1] 
  }) %>% unlist())

# plot at maximum 500 cells per group
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
ggsave("../../out/image/revision/119_ViolinCluster_technical_res0.1.pdf",width = 7,height = 12)

# pick one resolution onward 0.1
# main umap
plot03 <- DimPlot(data.combined, reduction = "umap", group.by = "RNA_snn_res.0.1",label = T,raster = T)
ggsave(plot=plot03,"../../out/image/revision/119_UMAPCluster_res0.1.pdf",width = 6,height = 5)

plot03a <- DimPlot(data.combined, reduction = "umap", group.by = "RNA_snn_res.0.1",label = T,raster = T,split.by = "orig.ident",ncol = 5)
ggsave(plot=plot03a,"../../out/image/revision/119_UMAPClusterSplit_res0.1.pdf",width = 13,height = 12)

# main umap cell cycle
plot04 <- DimPlot(data.combined, reduction = "umap", group.by = "Phase",raster = T,order = T)
ggsave(plot=plot04,"../../out/image/revision/119_UMAPPhase_res0.1.pdf",width = 5,height = 4)

# split by sample
Idents(data.combined) <- "RNA_snn_res.0.1"
df_meta <- data.combined@meta.data %>%
  rownames_to_column("rowname")
df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("rowname")

# plot the proporition for the phase per cluster
df_summary_phase <- df_meta %>%
  group_by(RNA_snn_res.0.1,Phase) %>%
  summarise(n = n()) %>%
  group_by(RNA_snn_res.0.1) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)

df_summary_phase %>%
  mutate(RNA_snn_res.0.1 = fct_relevel(RNA_snn_res.0.1,as.character(0:15))) %>%
  ggplot() +
  geom_col(aes(x=RNA_snn_res.0.1,y=prop,fill=Phase))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/image/revision/119_BarplotPhase_summary_res0.1.pdf",width = 7,height = 6)

# split by donor
df_summary_phase2 <- df_meta %>%
  group_by(orig.ident,RNA_snn_res.0.1,Phase) %>%
  summarise(n = n()) %>%
  group_by(orig.ident,RNA_snn_res.0.1) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)

df_summary_phase2 %>%
  mutate(RNA_snn_res.0.1 = fct_relevel(RNA_snn_res.0.1,as.character(0:15))) %>%
  ggplot() +
  geom_col(aes(x=orig.ident,y=prop,fill=Phase))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  facet_wrap(~RNA_snn_res.0.1,ncol = 1)+theme(strip.background = element_blank())
ggsave("../../out/image/revision/119_BarplotPhaseSplit_summary_res0.1.pdf",width = 9,height = 20)

# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(orig.ident,RNA_snn_res.0.1) %>%
  summarise(n = n()) %>%
  group_by(orig.ident) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary,"../../out/table/revision/119_summaryCluster_res0.1.tsv")

color_id <- alphabet(length(unique(df_summary$RNA_snn_res.0.1)))
# check the colors
show_col(color_id)

df_summary %>%
  ggplot() +
  geom_col(aes(x=orig.ident,y=prop,fill=RNA_snn_res.0.1))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  scale_fill_manual(values = unname(color_id))
ggsave("../../out/image/revision/119_Barplot_summary_res0.1.pdf",width = 10,height = 6)

# render the same plot as an heatmap
sample_prop_wide <- df_summary %>%
  # scale by rows
  group_by(RNA_snn_res.0.1) %>%
  mutate(zscore = (prop-mean(prop))/sd(prop)) %>%
  # make it long
  dplyr::select(orig.ident,RNA_snn_res.0.1,zscore) %>%
  pivot_wider(names_from = orig.ident,values_from = zscore,values_fill = 0) %>%
  column_to_rownames("RNA_snn_res.0.1")

rowSums(sample_prop_wide)
colSums(sample_prop_wide)

# plot the data as heatmap
meta_sample_prop <- data.frame(orig.ident = colnames(sample_prop_wide)) %>%
  left_join(df_meta %>%
              group_by(orig.ident,origin,sex,pathology_class) %>%
              summarise(),by=c("orig.ident"))

color_id2 <- alphabet(length(unique(meta_sample_prop$pathology_class)))
# check the colors
show_col(color_id2)

# build the named vector
names(color_id2) <- unique(meta_sample_prop$pathology_class)

column_meta_sample_prop <- HeatmapAnnotation(gender = meta_sample_prop$sex,
                                             origin = meta_sample_prop$origin,
                                             diagnosis = meta_sample_prop$pathology_class,
                                             col = list(gender = c("M" = "blue",
                                                                   "F" = "pink"),
                                                        origin = c("wm" = "gray",
                                                                   "cortex" = "black"),
                                                        diagnosis = color_id2))

ht2_shr_MG2 <- Heatmap(sample_prop_wide, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                       name = "zscore \nprop_cell_type \nscale cluster",
                       column_title = "sample",
                       # col = viridis::viridis(option = "turbo",n = 10),
                       
                       # row_names_gp = gpar(fontsize = 3),
                       top_annotation = column_meta_sample_prop,show_row_names = T
                       # cluster_rows = F,
                       # right_annotation = row_ha,
                       # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                       
)
pdf("../../out/image/revision/119_HeatmapCluster_summary_res0.1.pdf",width = 10,height = 6)
draw(ht2_shr_MG2,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# render the same plot as an heatmap
sample_prop_wide2 <- df_summary %>%
  # scale by rows
  group_by(orig.ident) %>%
  mutate(zscore = (prop-mean(prop))/sd(prop)) %>%
  # make it long
  dplyr::select(orig.ident,RNA_snn_res.0.1,zscore) %>%
  pivot_wider(names_from = orig.ident,values_from = zscore,values_fill = 0) %>%
  column_to_rownames("RNA_snn_res.0.1")

rowSums(sample_prop_wide2)
colSums(sample_prop_wide2)

meta_sample_prop2 <- data.frame(orig.ident = colnames(sample_prop_wide2)) %>%
  left_join(df_meta %>%
              group_by(orig.ident,origin,sex,pathology_class) %>%
              summarise(),by=c("orig.ident"))

# plot the data as heatmap
color_id3 <- alphabet(length(unique(meta_sample_prop2$pathology_class)))
# check the colors
show_col(color_id3)

# build the named vector
names(color_id3) <- unique(meta_sample_prop2$pathology_class)

column_meta_sample_prop2 <- HeatmapAnnotation(gender = meta_sample_prop2$sex,
                                              origin = meta_sample_prop2$origin,
                                              diagnosis = meta_sample_prop2$pathology_class,
                                             col = list(gender = c("M" = "blue",
                                                                   "F" = "pink"),
                                                        origin = c("wm" = "gray",
                                                                   "cortex" = "black"),
                                                        diagnosis = color_id3))

ht2_shr_MG22 <- Heatmap(sample_prop_wide2, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                        name = "zscore \nprop_cell_type \nscale sample",
                        column_title = "sample",
                        # col = viridis::viridis(option = "turbo",n = 10),
                        
                        # row_names_gp = gpar(fontsize = 3),
                        top_annotation = column_meta_sample_prop2,show_row_names = T
                        # cluster_rows = F,
                        # right_annotation = row_ha,
                        # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                        
)

pdf("../../out/image/revision/119_HeatmapCluster_summary2_res0.1.pdf",width = 10,height = 6)
draw(ht2_shr_MG22,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# # render the same plot as an heatmap
# sample_prop_wide3 <- df_summary_02 %>%
#   # make it long
#   dplyr::select(original_sample_name,RNA_snn_res.0.2,prop) %>%
#   pivot_wider(names_from = original_sample_name,values_from = prop,values_fill = 0) %>%
#   column_to_rownames("RNA_snn_res.0.2")
# 
# rowSums(sample_prop_wide3)
# colSums(sample_prop_wide3)
# 
# # plot the data as heatmap
# meta_sample_prop3 <- data.frame(colname = colnames(sample_prop_wide3)) %>%
#   left_join(df_meta %>%
#               group_by(original_sample_name,sex,diagnosis,location) %>%
#               summarise(),by=c("colname" = "original_sample_name"))
# 
# column_meta_sample_prop3 <- HeatmapAnnotation(gender = meta_sample_prop3$sex,
#                                               location = meta_sample_prop3$location,
#                                               diagnosis = meta_sample_prop3$diagnosis,
#                                               col = list(gender = c("m" = "blue",
#                                                                     "f" = "pink"),
#                                                          location = c("lumbar" = "gray90",
#                                                                       "thoracic" = "gray50",
#                                                                       "cervical" = "black"),
#                                                          diagnosis = c("Non-demented control" = "green",
#                                                                        "Multiple sclerosis" = "red")))
# 
# ht2_shr_MG23 <- Heatmap(sample_prop_wide3, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
#                         name = "prop_cell_type \nscale sample",
#                         column_title = "sample",
#                         col = viridis::viridis(option = "turbo",n = 10),
#                         
#                         # row_names_gp = gpar(fontsize = 3),
#                         top_annotation = column_meta_sample_prop3,show_row_names = T
#                         # cluster_rows = F,
#                         # right_annotation = row_ha,
#                         # row_split = rep(c(1,2,3,4),c(2,3,4,7))
#                         
# )
# pdf("../../out/plot/manualClean/HeatmapSampleProp_summary_harmonySkipIntegration_AllSoupX_00500_07000_05_res0.2.pdf",width = 10,height = 6)
# draw(ht2_shr_MG23,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
# dev.off()
# 
# # # define a convenient palette of colors
# # show_col(hue_pal()(9))
# # RColorBrewer::display.brewer.all()
# # col_pal <- RColorBrewer::brewer.pal(name = "Paired",n = 9)
# # # col_pal <- c("#E6E6E6","#ffff09","#c1ce08","#446d05","#053c03","#4D4D4D","#06a8ce","#033b6d","#ff0ed7","#9a0404")
# # show_col(col_pal)

# Identify conserved cell type markers ------------------------------------
# data
DefaultAssay(data.combined) <- "RNA"
# notice that in this case the data object of the RNA slot is already filled with the normalzied data, therefore in this case (following the Normalize workfloe for the integration) there is no need to run the NormalizeData on the RNA slot of the integrated object
# sobj_total_h@assays$RNA@data[20:50,1:10]
# dim(sobj_total_h@assays$RNA@data)
# 
# # scale the data see the note on evernote on why this can be done also at this point. this is needed because the scale.data is empty
# sobj_total_h@assays$RNA@scale.data
# all.genes <- rownames(sobj_total_h)
# sobj_total_h <- ScaleData(sobj_total_h,vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"))
# # confirm now the slot is loaded
# sobj_total_h@assays$RNA@scale.data[1:10,1:10]

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Idents(data.combined) <- "RNA_snn_res.0.1"
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/revision/119_FindAllMarkers_res0.1.tsv")

# save the top 100
sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  write_tsv("../../out/table/revision/119_FindAllMarkers_res0.1_top100.tsv")

# sobj_total_h.markers <- read_tsv("../../out/table/FindAllMarkers_harmonySkipIntegration_AllSoupX_01000_06000_15.tsv")

# try plotting the top markers
top_specific_markers <- sobj_total_h.markers %>%
  # filter ribosomal and mt genes
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC)

# And generate e.g. a dotplot:
dittoSeq::dittoDotPlot(data.combined,
                       vars = unique(top_specific_markers$gene), 
                       group.by = "RNA_snn_res.0.2")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
ggsave("../../out/image/revision/119_TopMarkersDitto_res0.1.pdf",width = 15,height = 6)

# # save the object ---------------------------------------------------------
# # save the object with the full annotation
# saveRDS(data.combined,"../../out/object/manualClean/data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_Annotation.rds")
# 
# # EDA ---------------------------------------------------------------------
# # compare some specific clusters
# Idents(data.combined) <- "RNA_snn_res.0.2"
# sobj_total_res0.2_9_0.markers <- RunPresto(data.combined, ident.1 = "9",ident.2 = "0",only.pos = F, min.pct = 0.25, logfc.threshold = 0.1)
# 
# # explore the HVG in the dataset
# VariableFeatures(data.combined) %>%
#   str_subset(pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
# 
# sobj_total_res0.2_9_0.markers %>%
#   ggplot(aes(x=avg_log2FC,y=-log(p_val_adj)))+geom_point()
# 
# library(enrichR)
# # DB selection ------------------------------------------------------------
# dbs <- listEnrichrDbs()
# # filter fo the db of interest
# dbs %>%
#   filter(str_detect(libraryName,pattern = "Atlas"))
# 
# dbs %>%
#   filter(str_detect(libraryName,pattern = "GO"))
# 
# dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023")
# 
# # query -------------------------------------------------------------------
# list_genes_UP <- list(clust_9_UP = sobj_total_res0.2_9_0.markers %>% 
#                         rownames_to_column("gene") %>%
#                         pull(gene))
# 
# # x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
# list_enrichr_UP <- lapply(list_genes_UP,function(x){
#   genes <- x
#   # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
#   out_enrich <- enrichr(genes, dbs_db)
#   
#   # filter out the annotations without an output
#   filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>%
#     unlist()
#   
#   out_enrich[filter_out>0] %>%
#     bind_rows(.id = "annotation")
# }) %>%
#   bind_rows(.id = "comparison")
# 
# list_enrichr_UP
# 
# 
# # list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")
# plot_list_UP <- list_enrichr_UP %>%
#   split(f = .$comparison)
# 
# # library(scales)
# list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
#   x %>%
#     group_by(annotation) %>%
#     arrange(P.value) %>%
#     dplyr::slice(1:10) %>%
#     mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
#     mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
#     # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
#     ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
#     scale_color_gradientn(colors = c("red","blue"),
#                           values = rescale(c(0,1)),
#                           limits = c(0,0.2))+
#     theme(strip.background = element_blank(),
#           panel.border = element_rect(colour = "black", fill = NA))+
#     ggtitle(y)
#   # scale_color_gradient(low = "red",high = "blue")
#   
#   #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
# })
# 
# wrap_plots(list_plot_UP,nrow = 1)
# # ggsave("../../out/image/enrichR_DE_pseudobulk_ASTRO_ctrl_refCX_UP_shr_filterExp.pdf",width = 6,height = 12,limitsize = FALSE)