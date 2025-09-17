# # AIM ---------------------------------------------------------------------
# # the aim is to plot the data after integration. in particular to plot the sognature score for Martina's cell assignament.
# 
# # libraries ---------------------------------------------------------------
# library(harmony)
# library(Seurat)
# library(dplyr)
# library(cowplot)
# library(tidyverse)
# library(ggrepel)
# library(scales)
# library(RColorBrewer)
# library(SeuratWrappers)
# library(dittoSeq)
# library(clustree)
# library(pals)
# library(patchwork)
# library(ComplexHeatmap)
# library(circlize)
# 
# # read in the data --------------------------------------------------------
# data.combined <- readRDS("../../out/object/129_MG_subcluster_HarmonySample.rds")
# # data.combined2 <- readRDS("../../out/object/100_MG_subcluster_HarmonyRun.rds")
# (DimPlot(data.combined,group.by = "dataset") + ggtitle("Harmony Sample"))
# (DimPlot(data.combined,group.by = "RNA_snn_res.0.9") + ggtitle("Harmony Sample"))
# 
# # load the signatures defined by Martina for the immune cell types
# # use the top 100 DEGs she has defined  (the file was provided by martina)
# df_LUT <- data.frame(type = c("homeostatic","MIMs foamy","MIMs iron","DC_2","DC_6","stressed MG"),
#                      cluster = c(0,1,8,2,6,3))
# 
# list_sig <- read_csv("../../data/Signature_IMM_subset_Absinta2019.csv") %>%
#   left_join(df_LUT,by = "cluster") %>%
#   mutate(type = case_when(is.na(type)~paste0("clu_",cluster),
#                           T~type)) %>%
#   split(.$type) %>%
#   lapply(function(x){
#     x %>%
#       pull(gene) %>%
#       unique()
#   })
# 
# # wrangling ---------------------------------------------------------------
# # use the new clustering provided by Martina
# meta_new <- data.combined@meta.data %>%
#   mutate(cluster_martina = case_when(RNA_snn_res.0.9 %in% c(18) ~ "plasmablast",
#                                      RNA_snn_res.0.9 %in% c(8) ~ "T-cells",
#                                      RNA_snn_res.0.9 %in% c(10) ~ "per-mac",
#                                      RNA_snn_res.0.9 %in% c(6,4) ~ "dendritic",
#                                      RNA_snn_res.0.9 %in% c(12,19,5) ~ "mye-MG",
#                                      RNA_snn_res.0.9 %in% c(1,11) ~ "cont",
#                                      RNA_snn_res.0.9 %in% c(7,16) ~ "MIMS-foamy",
#                                      RNA_snn_res.0.9 %in% c(9) ~ "MIMS-iron",
#                                      RNA_snn_res.0.9 %in% c(13) ~ "stress-MG",
#                                      RNA_snn_res.0.9 %in% c(0,3) ~ "homeo-MG",
#                                      RNA_snn_res.0.9 %in% c(15) ~ "cl9Absinta",
#                                      RNA_snn_res.0.9 %in% c(17) ~ "cl17",
#                                      RNA_snn_res.0.9 %in% c(14) ~ "cl14",
#                                      RNA_snn_res.0.9 %in% c(2) ~ "und"))
# 
# # count the numerosity of the cluster and reorder them
# meta_new_summary <- meta_new %>%
#   group_by(cluster_martina) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))
# 
# meta_new_final <- meta_new %>%
#   mutate(cluster_martina = factor(cluster_martina,levels = meta_new_summary$cluster_martina)) %>%
#   mutate(cluster_martina_num = as.numeric(cluster_martina)-1)
#   
# # add the new classification to the original folder
# data.combined$cluster_martina <- meta_new_final$cluster_martina
# data.combined$cluster_martina_num <- meta_new_final$cluster_martina_num
# (DimPlot(data.combined,group.by = "cluster_martina_num",label= T) + ggtitle("Harmony Sample"))
# 
# # save the plot without number and remove cluster 1, 10 and 12
# data.combined_subset_test <- subset(data.combined,subset = cluster_martina_num %in% c(0,2,3,4,5,6,7,8,9,11,13))
# DimPlot(data.combined_subset_test,group.by = "cluster_martina_num",label= F)
# ggsave("../../out/image/129_UMAP_SigAbsinta2019_IMMLYMsubset_highres_test.pdf",height = 5,width = 6)
# 
# # plots -------------------------------------------------------------------
# # score the signatures for senescence
# data.combined <- Seurat::AddModuleScore(data.combined,
#                                         features = list_sig,
#                                         name = "_score")
# 
# df_rename_long <- data.frame(names = data.combined@meta.data %>%
#                                colnames() %>%
#                                str_subset("_score"),
#                              rename = paste0("score_",names(list_sig)))
# 
# lookup_long <- df_rename_long$names
# names(lookup_long) <- df_rename_long$rename
# 
# # rename the columns
# data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup_long))
# 
# # plot the scores from AddModuleScore
# list_plot_02_long <- lapply(df_rename_long$rename,function(x){
#   plot <- FeaturePlot(data.combined,features = x,order = T,
#                       reduction = "umap",
#                       raster = T) + scale_color_viridis_c(option = "turbo")
#   return(plot)
# })
# 
# wrap_plots(list_plot_02_long)
# ggsave("../../out/image/129_UMAP_SigAbsinta2019_IMMLYMsubset_highres_split_test.pdf",width = 25,height = 20)
# 
# # same as above but as violin plot
# list_plot <- lapply(df_rename_long$rename, function(x){ 
#   test <- VlnPlot(object = data.combined,features = x, group.by = "cluster_martina_num",raster = T)
#   return(test)
# })
# 
# # make it a dataframe
# # x <- list_plot[[1]]
# df_violin <- lapply(list_plot,function(x){ 
#   df <- x[[1]]$data 
#   
#   # extract the name of the gene 
#   feature <- colnames(df)[1] 
#   
#   df %>% 
#     mutate(feature = feature) %>% 
#     setNames(c("value","ident","feature")) 
# }) %>% 
#   bind_rows()
# 
# head(df_violin) 
# 
# # # plot at maximum 1000 cells per group
# # set.seed(123)
# # df_plot_violin <- df_violin %>% 
# #   group_by(ident,feature) %>%
# #   sample_n(size = 150,replace = F) %>%
# #   ungroup()
# # 
# # df_plot_violin_summary <- df_plot_violin %>%
# #   group_by(feature) %>%
# #   summarise(med_score = median(value))
# # 
# # df_plot_violin %>%
# #   ggplot(aes(y = ident, x = value)) + 
# #   geom_violin(scale = "width")+ 
# #   #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) + 
# #   geom_point(position = position_jitter(width = 0.2),alpha = 0.05,size = 0.5) + 
# #   facet_wrap(~feature,nrow = 1,scales = "free") + 
# #   theme_bw() + 
# #   geom_vline(data = df_plot_violin_summary,aes(xintercept = med_score),linetype="dashed",col="red") +
# #   theme(strip.background = element_blank(),
# #         axis.text.x = element_text(hjust = 1,angle = 45))
# # ggsave("../../out/image/122_Violin_senescence_res0.2.pdf",width = 30,height = 4)
# 
# # plot the average score per signature per cluster as a heatmap
# df_senescence <- data.combined@meta.data %>%
#   dplyr::select(cluster_martina_num,contains("score_")) %>%
#   pivot_longer(names_to = "signature",values_to = "score",-cluster_martina_num) %>%
#   group_by(cluster_martina_num,signature) %>%
#   summarise(avg_score = mean(score),
#             med_score = median(score)) %>%
#   mutate(cluster_id = paste0("clu_",cluster_martina_num,"_custom")) %>%
#   ungroup()
# 
# mat_senescence_avg <- df_senescence %>%
#   group_by(signature) %>%
#   mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
#   ungroup() %>%
#   select(signature,cluster_id,scaled_avg_score) %>%
#   pivot_wider(names_from = cluster_id,values_from = scaled_avg_score) %>%
#   column_to_rownames("signature")
# 
# hm01 <- Heatmap(mat_senescence_avg,name = "scaled_avg_exp",col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")))
# 
# mat_senescence_avg %>%
#   rownames_to_column("signatures") %>%
#   write_tsv("../../out/table/129_tableScaled_SigAbisnta2019_res0.9_IMMLYMsubset_highres_test.tsv")
# 
# df_senescence %>%
#   write_tsv("../../out/table/129_tableUnScaled_SigAbisnta2019_res0.9_IMMLYMsubset_highres_test.tsv")
# 
# pdf("../../out/image/129_heatmap_SigAbisnta2019_res0.9_IMMLYMsubset_highres_test.pdf",width = 10,height = 5)
# draw(hm01,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
# dev.off()
# 
# # try also the scaling in the other direction, cluster-wise
# mat_senescence_avg2 <- df_senescence %>%
#   group_by(cluster_id) %>%
#   mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
#   ungroup() %>%
#   select(signature,cluster_id,scaled_avg_score) %>%
#   pivot_wider(names_from = signature,values_from = scaled_avg_score) %>%
#   column_to_rownames("cluster_id")
# 
# hm02 <- Heatmap(mat_senescence_avg2,name = "scaled_avg_exp")
# 
# pdf("../../out/image/129_heatmap_SigAbisnta2019_res0.9_IMMLYMsubset_highres2_test.pdf",width = 8,height = 6)
# draw(hm02,heatmap_legend_side = "left",padding = unit(c(4, 2, 2, 80), "mm"))
# dev.off()
# 
# # dotplot marker panel
# shortlist_features_list <- list(
#   IMMUNE = c("LYVE1", "CD163", "MRC1", "LINGO1", "HSPA1A", "MOBP", "CD83","HIF1A", "CD22", "VEGFA", "SOD1", "TREM2", "CX3CR1", "P2RY12","C3", "CSF1R", "CSF1R","CD74", "RUNX1", "C1QB", "PTPRC", "AIF1", "HLA-DRA", "TYROBP"),
#   B_CELLS = c("IGHG1", "CD38"),
#   T_CELLS =  c("SKAP1", "CD8A", "CD2")
# )
# 
# # plot the shortlisted feature per cluster
# # notice that this is done only on the subset of the young (control) cells
# test_long01 <- DotPlot(data.combined,
#                        features = unique(unlist(shortlist_features_list)),
#                        dot.scale = 8,
#                        cluster.idents = T,
#                        group.by = "cluster_martina_num") +
#   RotatedAxis() +
#   labs(title = "cluster_martina_num")+
#   theme(strip.text = element_text(angle = 90))
# 
# df_test <- lapply(shortlist_features_list,function(x){
#   test_long01$data %>% 
#     filter(features.plot %in% x)
# }) %>% 
#   bind_rows(.id = "cell_type")
# 
# # plot mapping radius
# df_test %>%
#   # force the order
#   # mutate(id = factor(id,levels = c(3,11,0,6,2,12,13,9,4,5,7,8,1,10))) %>% 
#   mutate(cell_type = factor(cell_type,levels = c("IMMUNE","B_CELLS","T_CELLS"))) %>% 
#   ggplot(aes(x = features.plot,y = id)) +
#   geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
#   scale_radius(range = c(0, 8)) +
#   facet_grid(~cell_type,scales = "free",space = "free")+
#   theme_cowplot()+
#   theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+
#   scale_color_gradient(low = "lightgrey",high = "blue")
# ggsave(plot=test_long01,"../../out/image/129_DotplotLong_custom_MG_subcluster_test.pdf",width = 10,height = 6)
