# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration. in particular to plot the sognature score for Martina's cell assignament.

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
library(circlize)

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/132_AST_subcluster_HarmonySample_martina.rds")
# data.combined2 <- readRDS("../../out/object/100_MG_subcluster_HarmonyRun.rds")
(DimPlot(data.combined,group.by = "dataset") + ggtitle("Harmony Sample"))

# load the signatures defined by Martina for the immune cell types
# use the top 100 DEGs she has defined  (the file was provided by martina)
# df_LUT <- data.frame(type = c("homeostatic","MIMs foamy","MIMs iron","DC_2","DC_6","stressed MG"),
#                      cluster = c(0,1,8,2,6,3))

list_sig <- read_csv("../../data/Siganture_AST_subset_Absinta2019.csv") %>%
  # left_join(df_LUT,by = "cluster") %>%
  # mutate(type = case_when(is.na(type)~paste0("clu_",cluster),
  #                         T~type)) %>%
  mutate(type = paste0("clu_Absinta2019_",cluster)) %>%
  split(.$type) %>%
  lapply(function(x){
    x %>%
      pull(gene) %>%
      unique()
  })

# plots -------------------------------------------------------------------
# score the signatures for senescence
data.combined <- Seurat::AddModuleScore(data.combined,
                                        features = list_sig,
                                        name = "_score")

df_rename_long <- data.frame(names = data.combined@meta.data %>%
                               colnames() %>%
                               str_subset("_score"),
                             rename = paste0("score_",names(list_sig)))

lookup_long <- df_rename_long$names
names(lookup_long) <- df_rename_long$rename

# rename the columns
data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup_long))

# plot the scores from AddModuleScore
list_plot_02_long <- lapply(df_rename_long$rename,function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

# wrap_plots(list_plot_02_long)
# ggsave("../../out/image/132_UMAP_SigAbsinta2019_ASTsubset.pdf",width = 25,height = 20)

# same as above but as violin plot
list_plot <- lapply(df_rename_long$rename, function(x){ 
  test <- VlnPlot(object = data.combined,features = x, group.by = "cellid",raster = T)
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

# # plot at maximum 1000 cells per group
# set.seed(123)
# df_plot_violin <- df_violin %>% 
#   group_by(ident,feature) %>%
#   sample_n(size = 150,replace = F) %>%
#   ungroup()
# 
# df_plot_violin_summary <- df_plot_violin %>%
#   group_by(feature) %>%
#   summarise(med_score = median(value))
# 
# df_plot_violin %>%
#   ggplot(aes(y = ident, x = value)) + 
#   geom_violin(scale = "width")+ 
#   #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) + 
#   geom_point(position = position_jitter(width = 0.2),alpha = 0.05,size = 0.5) + 
#   facet_wrap(~feature,nrow = 1,scales = "free") + 
#   theme_bw() + 
#   geom_vline(data = df_plot_violin_summary,aes(xintercept = med_score),linetype="dashed",col="red") +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/image/122_Violin_senescence_res0.2.pdf",width = 30,height = 4)

# plot the average score per signature per cluster as a heatmap
df_senescence <- data.combined@meta.data %>%
  dplyr::select(cellid,contains("score_")) %>%
  pivot_longer(names_to = "signature",values_to = "score",-cellid) %>%
  group_by(cellid,signature) %>%
  summarise(avg_score = mean(score),
            med_score = median(score)) %>%
  mutate(cluster_id = paste0(cellid)) %>%
  ungroup()

mat_senescence_avg <- df_senescence %>%
  group_by(signature) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = cluster_id,values_from = scaled_avg_score) %>%
  column_to_rownames("signature")

hm01 <- Heatmap(mat_senescence_avg,name = "scaled_avg_exp",col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

mat_senescence_avg %>%
  rownames_to_column("signatures") %>%
  write_tsv("../../out/table/132_tableScaled_SigAbisnta2019_cellid_ASTsubset.tsv")

df_senescence %>%
  write_tsv("../../out/table/132_tableUnScaled_SigAbisnta2019_cellid_ASTsubset.tsv")

pdf("../../out/image/132_heatmap_SigAbisnta2019_cellid_ASTsubset.pdf",width = 10,height = 5)
draw(hm01,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
dev.off()

# try also the scaling in the other direction, cluster-wise
mat_senescence_avg2 <- df_senescence %>%
  group_by(cluster_id) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = signature,values_from = scaled_avg_score) %>%
  column_to_rownames("cluster_id")

hm02 <- Heatmap(mat_senescence_avg2,name = "scaled_avg_exp")

pdf("../../out/image/132_heatmap_SigAbisnta2019_cellid_ASTsubset_2.pdf",width = 10,height = 5)
draw(hm02,heatmap_legend_side = "left",padding = unit(c(4, 2, 2, 80), "mm"))
dev.off()

# # keep only senamyo and plot the cluster score per disease area
# df_senescence2 <- data.combined@meta.data %>%
#   dplyr::select(RNA_snn_res.0.9,pathology_class,contains("senmayo")) %>%
#   group_by(RNA_snn_res.0.9,pathology_class) %>%
#   summarise(avg_score = mean(scoreSen_senmayo),
#             med_score = median(scoreSen_senmayo)) %>%
#   mutate(cluster_id = paste0("clu_",RNA_snn_res.0.9,"_res0.9")) %>%
#   ungroup()
# 
# mat_senescence_avg2 <- df_senescence2 %>%
#   group_by(cluster_id) %>%
#   mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
#   ungroup() %>%
#   select(pathology_class,cluster_id,scaled_avg_score) %>%
#   pivot_wider(names_from = pathology_class,values_from = scaled_avg_score) %>%
#   column_to_rownames("cluster_id")
# 
# hm02 <- Heatmap(mat_senescence_avg2,name = "scaled_avg_exp")
# 
# # pdf("../../out/image/129_heatmap_senescence_res0.9_IMMLYMsubset_highres.pdf",width = 12,height = 5)
# draw(hm02,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
# # dev.off()
# 
# # keep only senamyo and plot the cluster score per disease area
# mat_senescence_avg3 <- df_senescence2 %>%
#   group_by(pathology_class) %>%
#   mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
#   ungroup() %>%
#   select(pathology_class,cluster_id,scaled_avg_score) %>%
#   pivot_wider(names_from = cluster_id,values_from = scaled_avg_score) %>%
#   column_to_rownames("pathology_class")
# 
# hm03 <- Heatmap(mat_senescence_avg3,name = "scaled_avg_exp")
# 
# # pdf("../../out/image/129_heatmap_senescence_res0.9_IMMLYMsubset_highres.pdf",width = 12,height = 5)
# draw(hm03,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
# # dev.off()