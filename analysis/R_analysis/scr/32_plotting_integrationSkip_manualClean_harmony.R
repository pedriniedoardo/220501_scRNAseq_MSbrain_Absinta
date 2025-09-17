# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(GGally)
library(cowplot)

# read in the final object ------------------------------------------------
# scobj <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegration_AllSoupX.rds")
# this one is the dataset using 4000 HVG, which Martina recommended using.
scobj <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegration_AllSoupX_test.rds")

# is there the UMAP model in the object
scobj@reductions$umap@misc$model

# save the current meta add also the coordinates of the UMAP
df_umap <- scobj@reductions$umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column()

df_meta <- scobj@meta.data %>% 
  data.frame() %>% 
  rownames_to_column()

df_meta_full <- left_join(df_umap,df_meta,"rowname") %>% 
  mutate(orig_alt = case_when(orig.ident %in% c("s31","s27","s26","s25")~"wm_new",
                              T ~ origin))

df_meta_full %>%
  write_tsv("../../out/table/ManualClean/meta_data.combined_WM_CX_harmonySkipIntegAllSoupX.tsv")

# add also the varaible to the object
scobj$origin_alt <- df_meta_full$orig_alt

# check that all the samples are correctly assigned
# df_meta_full %>% 
#   group_by(origin,orig.ident) %>% 
#   summarise() %>% 
#   print(n=50)

df_meta_full %>% 
  group_by(origin,orig.ident) %>% 
  summarise() %>% 
  summarise(n = n())

# plotting ----------------------------------------------------------------
# confirm the identity of the object
DimPlot(scobj,label = T,raster = T)
ggsave("../../out/image/ManualClean/UMAP_data.combined_WM_CX_harmonySkipIntegAllSoupX.pdf",width = 8,height = 7)

DimPlot(scobj,label = T,raster = T,split.by = "orig.ident",ncol = 7)
ggsave("../../out/image/ManualClean/UMAP_data.combined_WM_CX_harmonySkipIntegAllSoupX_splitSample.pdf",width = 29,height = 28)

DimPlot(scobj,label = T,raster = T,split.by = "origin_alt")
ggsave("../../out/image/ManualClean/UMAP_data.combined_WM_CX_harmonySkipIntegAllSoupX_origin.pdf",width = 22,height = 7)

# add the splitting for the regular origin
DimPlot(scobj,label = T,raster = T,split.by = "origin")
ggsave("../../out/image/ManualClean/UMAP_data.combined_WM_CX_harmonySkipIntegAllSoupX_origin_regular.pdf",width = 15,height = 7)

DimPlot(scobj,label = T,raster = T,split.by = "origin_alt",group.by = "seurat_clusters_ref")
ggsave("../../out/image/ManualClean/UMAP_data.combined_WM_CX_harmonySkipIntegAllSoupX_originalMarkersID.pdf",width = 22,height = 7)

# df_meta_full %>% 
#   filter(origin=="wm") %>% 
#   ggplot(aes(x=UMAP_1,y=UMAP_2))+geom_point(size = 0.1)+
#   facet_wrap(~orig.ident)+theme_bw()+theme(strip.background = element_blank())
# 
# df_meta_full %>% 
#   filter(origin=="wm") %>% 
#   ggplot(aes(x=UMAP_1,y=UMAP_2))+geom_point(size = 0.1)+
#   facet_wrap(~seurat_clusters)+theme_bw()+theme(strip.background = element_blank())
# 
# df_meta_full %>% 
#   filter(origin=="wm") %>% 
#   ggplot(aes(x=UMAP_1,y=UMAP_2))+geom_point(size = 0.1)+
#   facet_grid(orig.ident~seurat_clusters)+theme_bw()+theme(strip.background = element_blank())
# 
# df_meta_full %>% 
#   filter(origin=="wm") %>% 
#   ggplot(aes(x=UMAP_1,y=UMAP_2))+geom_point(size = 0.1)+
#   facet_wrap(~seurat_clusters_ref)+theme_bw()+theme(strip.background = element_blank())

df_meta_full %>% 
  filter(origin=="wm") %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2))+geom_point(size = 0.1,alpha=0.5)+
  facet_grid(paste0("new_",seurat_clusters)~paste0("old_",seurat_clusters_ref))+theme_bw()+theme(strip.background = element_blank())

df_meta_full %>% 
  filter(origin=="cortex") %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2))+geom_point(size = 0.1,alpha=0.5)+
  facet_grid(paste0("new_",seurat_clusters)~paste0("old_",seurat_clusters_ref))+theme_bw()+theme(strip.background = element_blank())

# plot some metadata
# df_meta_full %>%
#   arrange(percent.mt.harmony) %>% 
#   ggplot(aes(x=UMAP_1,y=UMAP_2,col=percent.mt.harmony))+geom_point(size=0.1)+
#   scale_color_viridis_c(option = "turbo")+
#   theme_cowplot()
# ggsave("../../out/image/ManualClean/UMAP_data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_MTPercent.pdf",width = 8,height = 7)

df_meta_full %>%
  arrange(percent.mt) %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=percent.mt))+geom_point(size=0.1)+
  scale_color_viridis_c(option = "turbo")+
  theme_cowplot()

# ggsave("../../out/image/ManualClean/UMAP_data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_MTPercent.png",width = 8,height = 7)

# df_meta_full %>%
#   arrange(percent.mt.harmony) %>% 
#   ggplot(aes(x=UMAP_1,y=UMAP_2))+
#   stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE,n = 200) +
#   theme_cowplot()+
#   scale_fill_viridis_c(option = "turbo")+
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0))
# ggsave("../../out/image/ManualClean/UMAP_data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_MTPercentHM.pdf",width = 8,height = 7)

# save the meta for the proportion
df_summary <- scobj@meta.data %>% 
  group_by(orig.ident,origin,disease,pathology_class,seurat_clusters) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by() %>% 
  group_by(orig.ident) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop = n/tot)
write_tsv(df_summary,"../../out/table/ManualClean/df_summary_data.combined_WM_CX_harmonySkipIntegAllSoupX.tsv")

# plot the propotions per cluster
df_summary %>% 
  ggplot(aes(x=factor(seurat_clusters),y=prop,col=disease)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.6,jitter.width = 0.1))+theme_bw()+facet_wrap(~origin,nrow = 2)+theme(strip.background = element_blank())+scale_y_log10()
ggsave("../../out/image/ManualClean/plot_summary_data.combined_WM_CX_harmonySkipIntegAllSoupX_log.pdf",width = 8,height = 8)

df_summary %>% 
  ggplot(aes(x=factor(seurat_clusters),y=prop,col=disease)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.6,jitter.width = 0.1))+theme_bw()+facet_wrap(~origin,nrow = 2)+theme(strip.background = element_blank())
ggsave("../../out/image/ManualClean/plot_summary_data.combined_WM_CX_harmonySkipIntegAllSoupX_lin.pdf",width = 8,height = 8)

# do the same for the proportion per tissue
df_summary %>% 
  ggplot(aes(x=origin,y=prop,col=disease)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.6,jitter.width = 0.1))+theme_bw()+facet_wrap(~seurat_clusters,scales = "free")+theme(strip.background = element_blank())
ggsave("../../out/image/ManualClean/plot_summary_data.combined_WM_CX_harmonySkipIntegAllSoupX_log.pdf",width = 8,height = 8)

# test <- scobj@meta.data %>% 
#   data.frame() %>% 
#   dplyr::select(seurat_clusters,nCount_RNA.harmony,nFeature_RNA.harmony,Phase,percent.mt.harmony,percent.ribo)
# 
# # plot all the metadata
# test %>% 
#   rownames_to_column("barcode") %>% 
#   pivot_longer(names_to="variable",values_to="value",-c(seurat_clusters,barcode,Phase)) %>% 
#   ggplot(aes(x=seurat_clusters,y=value))+geom_violin()+geom_boxplot(width=0.1,outlier.shape = NA)+facet_wrap(~variable,scales = "free")+theme_bw()+theme(strip.background = element_blank())+ scale_y_log10()
# ggsave("../../out/image/ManualClean/plot_metrics_manualClean_SoupX_harmony.pdf",width = 8,height = 6)

# plot all the 
# ggpairs(test[,-1],)

# save the panel martina uses for the annotation
shortlist_features_list_long <- list(
  IMMUNE = c("IGHG1","CD38","CD8A","CD2","SKAP1","LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
  NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
  PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1")
)

shortlist_features_martina <- list(
  LAYERS = c("ADRA2A", "NR4A2", "NPY", "NOS1", "SST", "SLC6A1", "SLC6A12","NDNF","RASGRF2","CUX2","RORB","TRIB2","B3GALT2","NTNG2","TLE4","CCN2","CALB2","ERBB4","PVALB","TAC3","VIP")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long2 <- DotPlot(scobj, features = shortlist_features_list_long, dot.scale = 8,cluster.idents = T)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# theme(axis.text.x = element_text(angle = 45,hjust = 1))
test_long2
ggsave("../../out/image/ManualClean/Dotplot_data.combined_WM_CX_harmonySkipIntegAllSoupX.pdf",width = 23,height = 4)

# test_martina <- DotPlot(scobj, features = shortlist_features_martina, dot.scale = 8,cluster.idents = T) +
#   RotatedAxis()
# test_martina
# ggsave("out/image/Dotplot_clusters_NOT_annotated_big_panle_integration_norm_fix_martina_regressCC_DoubletSinglet.pdf",width = 10,height = 4)

# -------------------------------------------------------------------------
# show the original UMAPs from the individual objects
CX_ansinta <- readRDS("../../out/object/ManualClean/data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony.rds")
DimPlot(CX_ansinta,raster=T,label=T)

WM_ansinta <- readRDS("../../data/WM_data/WM_data_Absinta_Nature2021/all20_integrated_clean_metadata.rds")
DimPlot(WM_ansinta,raster=T,label=T,reduction = "curated")


# tailored plotting -------------------------------------------------------
# Hi edo, can you  print a simplified version of the dotplot?
# IMMUNE: CX3CR1, P2RY12, CSF1R, RUNX1, SKAP1, PTPRC
# OLIGO: PLP1, MBP, MOG
# OPC: SOX6, PDGFRA, NLGN4X
# ASTRO: AQP4, SCL1A2, GFAP
# VAS: CLDN5, FLT1
# NEURONS: SYT1, GAD2, RORB, SATB2, CUX2, NRGN, SCL17A7
# Ordina cluster cosi’
# 5, 13, 0, 4, 14, 9, 3, 12, 11, 7, 8, 10, 6, 2, 1
shortlist_features3 <- list(
  IMMUNE = c("CX3CR1", "P2RY12", "CSF1R", "RUNX1", "SKAP1", "PTPRC"),
  OLIGO = c("PLP1", "MBP", "MOG"),
  OPC = c("SOX6", "PDGFRA", "NLGN4X"),
  ASTRO = c("AQP4", "SLC1A2", "GFAP"),
  VAS = c("CLDN5", "FLT1"),
  NEURONS = c("SYT1", "GAD2", "RORB", "SATB2", "CUX2", "NRGN", "SLC17A7")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long3 <- DotPlot(scobj, features = shortlist_features3, dot.scale = 8,cluster.idents = T)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# test_long3
# ggsave("../../out/image/ManualClean/Dotplot_data.combined_WM_CX_harmonySkipIntegAllSoupX_tailored.pdf",width = 12,height = 4)
# force the order suggested by matrina
df_test <- lapply(shortlist_features3,function(x){
  test_long3$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

glimpse(df_test)

df_test %>%
  # force the order
  mutate(id = factor(id,levels = c(5, 13, 0, 4, 14, 9, 3, 12, 11, 7, 8, 10, 6, 2, 1))) %>% 
  mutate(cell_type = factor(cell_type,levels = c("IMMUNE","OLIGO","OPC","ASTRO","VAS","NEURONS"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_size(range = c(0, 6)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")
ggsave("../../out/image/ManualClean/Dotplot_data.combined_WM_CX_harmonySkipIntegAllSoupX_tailoredReorder.pdf",width = 8,height = 4)

# for the file with the proportion, can you please provide a new one merging the following clusters?
# OLIGO: 0,4,14
# OPC: 9
# ASTRO: 3, 12
# IMMUNE: 5
# LYM: 13
# VAS: 11
# EXC NEU: 1, 2, 10,6
# INH NEU: 7, 8
df_summary2 <- scobj@meta.data %>% 
  mutate(cell_id = case_when(seurat_clusters %in% c(0,4,14)~"OLIGO",
                             seurat_clusters %in% c(9)~"OPC",
                             seurat_clusters %in% c(3,12)~"ASTRO",
                             seurat_clusters %in% c(5)~"IMMUNE",
                             seurat_clusters %in% c(13)~"LYM",
                             seurat_clusters %in% c(11)~"VAS",
                             seurat_clusters %in% c(1, 2, 10,6)~"EXC NEU",
                             seurat_clusters %in% c(7,8)~"INH NEU")) %>% 
  group_by(orig.ident,origin,disease,pathology_class,cell_id) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by() %>% 
  group_by(orig.ident) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop = n/tot)
write_tsv(df_summary2,"../../out/table/ManualClean/df_summary_data.combined_WM_CX_harmonySkipIntegAllSoupX_cellID.tsv")

df_summary2 %>% 
  ggplot(aes(x=origin,y=prop,col=disease)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.6,jitter.width = 0.1))+
  theme_bw()+
  facet_wrap(~cell_id,nrow = 2,scales = "free")+
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("blue","red"))
ggsave("../../out/image/ManualClean/plot_summary_data.combined_WM_CX_harmonySkipIntegAllSoupX_CellID.pdf",width = 12,height = 8)
