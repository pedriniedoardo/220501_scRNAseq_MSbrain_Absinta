# # libraries ---------------------------------------------------------------
# library(tidyverse)
# library(speckle)
# library(limma)
# library(statmod)
# library(cowplot)
# library(ggrepel)
# 
# # # read in the data --------------------------------------------------------
# # so <- readRDS("../out_large/scRNAseq_analysis/object/sobj_total_h_fix_filter_norm_doublet_harmony_5K_scaleSbatch_SCtypeAnnotation.rds")
# # DimPlot(so,label = T,group.by = "annotation_confident")
# # 
# # # load the annotation for the subclusters
# # meta_subNeu <- read_tsv("out/table/meta_subNeu_refHippo_classification.tsv") %>% 
# #   select(barcodes,contains("assigned"))
# # 
# # # wrangling ---------------------------------------------------------------
# # # add the annotation to the dataset
# # meta_ref <- so@meta.data %>% 
# #   rownames_to_column("barcodes") %>% 
# #   left_join(meta_subNeu,by = "barcodes") %>% 
# #   mutate(assigned.predicted.id.l1 = case_when(is.na(assigned.predicted.id.l1)~annotation_confident,
# #                                               T~assigned.predicted.id.l1),
# #          assigned.predicted.id.l2 = case_when(is.na(assigned.predicted.id.l2)~annotation_confident,
# #                                               T~assigned.predicted.id.l2),
# #          assigned.predicted.id.l3 = case_when(is.na(assigned.predicted.id.l3)~annotation_confident,
# #                                               T~assigned.predicted.id.l3)) %>% 
# #   column_to_rownames("barcodes")
# # 
# # # update the meta
# # so@meta.data <- meta_ref
# # 
# # # save the new heatmap
# # DimPlot(so,group.by = "assigned.predicted.id.l1",label = T)
# # 
# # # save the meta
# # write_tsv(meta_ref %>% rownames_to_column("barcodes"),file = "out/table/meta_sobj_total_h_fix_filter_norm_doublet_harmony_5K_scaleSbatch_SCtypeAnnotation_HippoClass.tsv")
# 
# meta_ref <- read_tsv(file = "out/table/meta_sobj_total_h_fix_filter_norm_doublet_harmony_5K_scaleSbatch_SCtypeAnnotation_HippoClass.tsv")
# 
# # explore the data --------------------------------------------------------
# # confirm the numbers from the tissue dataset
# meta_ref %>% 
#   group_by(ADNPscoring,assigned.predicted.id.l1) %>% 
#   summarise(n = n())
# 
# # brain MS vs control -----------------------------------------------------
# # run it on the predictred l2, not the robust one to avoid the missing values
# # Run propeller testing for cell type proportion differences between the groups
# out_brain <- propeller(clusters = meta_ref$assigned.predicted.id.l1,
#                        sample = meta_ref$orig.ident.cca,
#                        group = meta_ref$ADNPscoring)
# 
# out_brain %>% 
#   rownames_to_column("cluster") %>% 
#   write_tsv("out/table/propeller_out.tsv")
# 
# 
# # plotting
# # plotCellTypeProps(clusters = meta_ref$assigned.predicted.id.l1,
# #                   sample = meta_ref$orig.ident.cca)
# df_summary <- meta_ref %>% 
#   group_by(assigned.predicted.id.l1,
#            orig.ident.cca,
#            ADNPscoring) %>% 
#   summarise(n=n()) %>% 
#   ungroup() %>% 
#   group_by(orig.ident.cca) %>% 
#   mutate(tot = sum(n),
#          prop = n/tot)
# 
# # plot
# df_summary %>% 
#   ggplot(aes(x=ADNPscoring,y=prop))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1))+
#   facet_wrap(~assigned.predicted.id.l1,scales = "free")+
#   theme_bw()+
#   theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("out/image/propeller_plot.pdf",width = 20,height = 15)
# 
# df_summary %>% 
#   ggplot(aes(x=ADNPscoring,y=prop))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1))+
#   geom_text_repel(aes(label = orig.ident.cca,x=ADNPscoring,y=prop))+
#   facet_wrap(~assigned.predicted.id.l1,scales = "free")+
#   theme_bw()+
#   theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("out/image/propeller_plot2.pdf",width = 20,height = 15)
# 
# df_summary %>%
#   ungroup() %>% 
#   group_by(assigned.predicted.id.l1,ADNPscoring) %>% 
#   summarise(n = sum(n)) %>% 
#   ungroup() %>% 
#   group_by(ADNPscoring) %>% 
#   mutate(tot = sum(n),
#          prop = n/tot) %>% 
#   ggplot(aes(x=ADNPscoring,y=prop,fill=assigned.predicted.id.l1))+
#   geom_bar(stat = "identity", position="fill") +
#   ggtitle("SCType predictions") +
#   theme_cowplot() +
#   labs(x = "SampleID", y = "fraction of cells") 
# ggsave("out/image/propeller_plot_stacked.pdf",width = 8,height = 6)
# 
# 
# # -------------------------------------------------------------------------
# # test remove sample 10
# meta_ref <- read_tsv(file = "out/table/meta_sobj_total_h_fix_filter_norm_doublet_harmony_5K_scaleSbatch_SCtypeAnnotation_HippoClass.tsv") %>% 
#   filter(orig.ident.cca != "10_cr_61")
# 
# # explore the data --------------------------------------------------------
# # confirm the numbers from the tissue dataset
# meta_ref %>% 
#   group_by(ADNPscoring,assigned.predicted.id.l1) %>% 
#   summarise(n = n())
# 
# # brain MS vs control -----------------------------------------------------
# # run it on the predictred l2, not the robust one to avoid the missing values
# # Run propeller testing for cell type proportion differences between the groups
# out_brain <- propeller(clusters = meta_ref$assigned.predicted.id.l1,
#                        sample = meta_ref$orig.ident.cca,
#                        group = meta_ref$ADNPscoring)
# 
# out_brain %>% 
#   rownames_to_column("cluster") %>% 
#   write_tsv("out/table/propeller_out_remove10.tsv")
# 
# 
# # plotting
# # plotCellTypeProps(clusters = meta_ref$assigned.predicted.id.l1,
# #                   sample = meta_ref$orig.ident.cca)
# df_summary <- meta_ref %>% 
#   group_by(assigned.predicted.id.l1,
#            orig.ident.cca,
#            ADNPscoring) %>% 
#   summarise(n=n()) %>% 
#   ungroup() %>% 
#   group_by(orig.ident.cca) %>% 
#   mutate(tot = sum(n),
#          prop = n/tot)
# 
# df_summary %>% 
#   ggplot(aes(x=ADNPscoring,y=prop))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1))+
#   geom_text_repel(aes(label = orig.ident.cca,x=ADNPscoring,y=prop))+
#   facet_wrap(~assigned.predicted.id.l1,scales = "free")+
#   theme_bw()+
#   theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("out/image/propeller_plot2_remove10.pdf",width = 20,height = 15)
# 
# # -------------------------------------------------------------------------
# # test remove sample 3
# meta_ref <- read_tsv(file = "out/table/meta_sobj_total_h_fix_filter_norm_doublet_harmony_5K_scaleSbatch_SCtypeAnnotation_HippoClass.tsv") %>% 
#   filter(orig.ident.cca != "03_cr_61")
# 
# # explore the data --------------------------------------------------------
# # confirm the numbers from the tissue dataset
# meta_ref %>% 
#   group_by(ADNPscoring,assigned.predicted.id.l1) %>% 
#   summarise(n = n())
# 
# # brain MS vs control -----------------------------------------------------
# # run it on the predictred l2, not the robust one to avoid the missing values
# # Run propeller testing for cell type proportion differences between the groups
# out_brain <- propeller(clusters = meta_ref$assigned.predicted.id.l1,
#                        sample = meta_ref$orig.ident.cca,
#                        group = meta_ref$ADNPscoring)
# 
# out_brain %>% 
#   rownames_to_column("cluster") %>% 
#   write_tsv("out/table/propeller_out_remove03.tsv")
# 
# 
# # plotting
# # plotCellTypeProps(clusters = meta_ref$assigned.predicted.id.l1,
# #                   sample = meta_ref$orig.ident.cca)
# df_summary <- meta_ref %>% 
#   group_by(assigned.predicted.id.l1,
#            orig.ident.cca,
#            ADNPscoring) %>% 
#   summarise(n=n()) %>% 
#   ungroup() %>% 
#   group_by(orig.ident.cca) %>% 
#   mutate(tot = sum(n),
#          prop = n/tot)
# 
# df_summary %>% 
#   ggplot(aes(x=ADNPscoring,y=prop))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1))+
#   geom_text_repel(aes(label = orig.ident.cca,x=ADNPscoring,y=prop))+
#   facet_wrap(~assigned.predicted.id.l1,scales = "free")+
#   theme_bw()+
#   theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("out/image/propeller_plot2_remove03.pdf",width = 20,height = 15)
