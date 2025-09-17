# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(GGally)
library(cowplot)

# read in the final object ------------------------------------------------
scobj <- readRDS("../../out/object/ManualClean/data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony.rds")

# save the current meta add also the coordinates of the UMAP
df_umap <- scobj@reductions$umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column()

df_meta <- scobj@meta.data %>% 
  data.frame() %>% 
  rownames_to_column()

df_meta_full <- left_join(df_umap,df_meta,"rowname")

df_meta_full %>%
  write_tsv("../../out/table/ManualClean/meta_data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony.tsv")

# plotting ----------------------------------------------------------------
# confirm the identity of the object
DimPlot(scobj,label = T,raster = T)
ggsave("../../out/image/ManualClean/UMAP_data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony.pdf",width = 8,height = 7)

DimPlot(scobj,label = T,raster = T,split.by = "official_id",ncol = 5)
ggsave("../../out/image/ManualClean/UMAP_data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_splitSample.pdf",width = 26,height = 20)

# plot some metadata
# df_meta_full %>%
#   arrange(percent.mt.harmony) %>% 
#   ggplot(aes(x=UMAP_1,y=UMAP_2,col=percent.mt.harmony))+geom_point(size=0.1)+
#   scale_color_viridis_c(option = "turbo")+
#   theme_cowplot()
# ggsave("../../out/image/ManualClean/UMAP_data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_MTPercent.pdf",width = 8,height = 7)

df_meta_full %>%
  arrange(percent.mt.harmony) %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=percent.mt.harmony))+geom_point(size=0.1)+
  scale_color_viridis_c(option = "turbo")+
  theme_cowplot()
ggsave("../../out/image/ManualClean/UMAP_data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony_MTPercent.png",width = 8,height = 7)

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
  group_by(official_id,orig.ident.cca,disease,seurat_clusters) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by() %>% 
  group_by(official_id,orig.ident.cca,disease) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop = n/tot)
write_tsv(df_summary,"../../out/table/ManualClean/df_summary_manualClean_SoupX_harmony.tsv")

# plot the propotions per cluster
df_summary %>% 
  ggplot(aes(x=seurat_clusters,y=prop,col=disease)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.6,jitter.width = 0.1))+theme_bw()
ggsave("../../out/image/ManualClean/plot_summary_manualClean_SoupX_harmony.pdf",width = 8,height = 4)

test <- scobj@meta.data %>% 
  data.frame() %>% 
  dplyr::select(seurat_clusters,nCount_RNA.harmony,nFeature_RNA.harmony,Phase.harmony,percent.mt.harmony,percent.ribo.harmony)

# plot all the metadata
test %>% 
  rownames_to_column("barcode") %>% 
  pivot_longer(names_to="variable",values_to="value",-c(seurat_clusters,barcode,Phase.harmony)) %>% 
  ggplot(aes(x=seurat_clusters,y=value))+geom_violin()+geom_boxplot(width=0.1,outlier.shape = NA)+facet_wrap(~variable,scales = "free")+theme_bw()+theme(strip.background = element_blank())+ scale_y_log10()
ggsave("../../out/image/ManualClean/plot_metrics_manualClean_SoupX_harmony.pdf",width = 8,height = 6)

# plot all the 

ggpairs(test[,-1],)

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
ggsave("../../out/image/ManualClean/Dotplot_data.combined_fix_filter_norm_doublet_integrated_manualClean_SoupX_harmony.pdf",width = 23,height = 4)

# test_martina <- DotPlot(scobj, features = shortlist_features_martina, dot.scale = 8,cluster.idents = T) +
#   RotatedAxis()
# test_martina
# ggsave("out/image/Dotplot_clusters_NOT_annotated_big_panle_integration_norm_fix_martina_regressCC_DoubletSinglet.pdf",width = 10,height = 4)
