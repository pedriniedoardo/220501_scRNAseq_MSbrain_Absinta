# AIM ---------------------------------------------------------------------
# load the final object and print the tables needed for martina


# libraries ---------------------------------------------------------------
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
scobj <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")

# summary table per sample ------------------------------------------------
# generate a table of numbe rof nuclei per saample per annotation
df_summary1 <- scobj@meta.data %>%
  group_by(orig.ident,sample_id,facility,origin,disease,sex,age,pathology,pathology_class,patient,plaque,PMI) %>%
  summarise(mean_nCount = mean(nCount_RNA),
            mean_nFeature = mean(nFeature_RNA),
            n_nuclei = n())

# make a summary also for the number of cells based on the annotation
df_summary2 <- scobj@meta.data %>%
  group_by(orig.ident,expertAnno.l1) %>%
  summarise(n_cell = n()) %>%
  pivot_wider(names_from = expertAnno.l1,values_from = n_cell,values_fill = 0)

# join the tables
df_summary_tot <- df_summary1 %>%
  left_join(df_summary2,by = "orig.ident")

write_tsv(df_summary_tot,file = "../../out/table/revision/120_extended_table_2.tsv")

# identify the top markers using the manual annotation --------------------
# data
DefaultAssay(scobj) <- "RNA"

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Idents(scobj) <- "expertAnno.l1"
sobj_total_h.markers <- RunPrestoAll(scobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/revision/120_FindAllMarkers_expertAnno.l1.tsv")

# save the top 100
sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  write_tsv("../../out/table/revision/120_FindAllMarkers_expertAnno.l1_top100.tsv")

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
                       group.by = "expertAnno.l1")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
ggsave("../../out/image/revision/120_TopMarkersDitto_expertAnno.l1.pdf",width = 15,height = 6)

# summary per pathologiocal stage -----------------------------------------
# summary table for the number of cells per pathological stage
# plot the proporition for the phase per cluster
df_summary_pahtology <- scobj@meta.data %>%
  group_by(pathology_class,expertAnno.l1) %>%
  summarise(n = n()) %>%
  group_by(pathology_class) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)

# save the table
df_summary_pahtology %>%
  write_tsv("../../out/table/revision/120_df_summary_pahtology.tsv")

df_summary_pahtology %>%
  ggplot() +
  geom_col(aes(x=pathology_class,y=prop,fill=expertAnno.l1))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/image/revision/120_Barplotpathology_summary_expertAnno.l1.pdf",width = 7,height = 6)


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

