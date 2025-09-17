# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration of the subclusters

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
data.combined <- readRDS("../../out/object/129_MG_subcluster_HarmonySample.rds")
# data.combined2 <- readRDS("../../out/object/100_MG_subcluster_HarmonyRun.rds")

(DimPlot(data.combined,group.by = "dataset") + ggtitle("Harmony Sample"))
(DimPlot(data.combined,group.by = "dataset",split.by = "origin") + ggtitle("Harmony Sample"))

# plots -------------------------------------------------------------------
# plot the tree of the cluster dependencies. this will justify the choice of the resolution, not too granula for the moment.
# library(clustree)
# clustree::clustree(data.combined@meta.data[,grep("RNA_snn_res", colnames(data.combined@meta.data))],
#                    prefix = "RNA_snn_res.")
# ggsave("../../out/image/129_UMAPCluster_tree_MG_subcluster.pdf",width = 10,height = 10)

# general UMAP with new clustering
DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.9",label = T,raster = T)
ggsave("../../out/image/129_UMAPCluster_res0.9_MG_subcluster.pdf",width = 6,height = 5)

DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.9",split.by = "origin",label = T,raster = T,ncol=4)
ggsave("../../out/image/129_UMAPCluster_splitTreat_res0.9_MG_subcluster.pdf",width = 8,height = 5)
# DimPlot(data.combined, reduction = "umap", split.by = "origin",label = T,raster = T,ncol=5)
# general UMAP with former clustering

# martina asked to split by both area and disease
df_all <- data.combined@meta.data %>%
  data.frame() %>%
  rownames_to_column() %>%
  # filter(! RNA_snn_res.0.9 %in% c(3)) %>%
  left_join(data.combined@reductions$umap@cell.embeddings %>%
              data.frame() %>%
              rownames_to_column(),by = "rowname")

# how many celle per condition
df_all %>%
  group_by(origin,disease) %>%
  summarise(n = n())

# subset the cells per area condition
obj_crop <- subset(data.combined, subset = RNA_snn_res.0.9 %in% c(0,2,3,9,7,16,13,4,6,10))
df_all_crop <- obj_crop@meta.data %>%
  data.frame() %>%
  rownames_to_column() %>%
  # filter(! RNA_snn_res.0.9 %in% c(3)) %>%
  left_join(data.combined@reductions$umap@cell.embeddings %>%
              data.frame() %>%
              rownames_to_column(),by = "rowname")
df_all_crop %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2)) +
  geom_point(size = 0.1,aes(col=RNA_snn_res.0.9)) +
  # geom_point(data = df_cluster,size = 0.1,col="#60b1bd") +
  theme_cowplot() +
  facet_grid(origin~disease)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(strip.background = element_blank())

# df_cluster <- data.combined@meta.data %>%
#   data.frame() %>%
#   rownames_to_column() %>%
#   filter(RNA_snn_res.0.9 %in% c(3)) %>%
#   left_join(data.combined@reductions$umap@cell.embeddings %>%
#               data.frame() %>%
#               rownames_to_column(),by = "rowname")
# 

df_all %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2)) +
  geom_point(size = 0.1,aes(col=RNA_snn_res.0.9)) +
  # geom_point(data = df_cluster,size = 0.1,col="#60b1bd") +
  theme_cowplot() +
  facet_grid(origin~disease)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(strip.background = element_blank())
ggsave("../../out/image/129_UMAPCluster_splitTreatOrigin_res0.9_MG_subcluster.pdf",width = 9,height = 8)

df_all %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2)) +
  geom_point(size = 0.1,aes(col=RNA_snn_res.0.9)) +
  # geom_point(data = df_cluster,size = 0.1,col="#60b1bd") +
  theme_cowplot() +
  facet_grid(origin~disease)+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(strip.background = element_blank())

# 
shortlist_features_list_long <- list(
  IMMUNE = c("LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC","AIF1","HLA-DRA","TYROBP"),
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
                       group.by = "RNA_snn_res.0.9") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.9")+
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_long01,"../../out/image/129_DotplotLong_res0.9_MG_subcluster.pdf",width = 30,height = 6)

# marker per cluster ------------------------------------------------------
DefaultAssay(data.combined) <- "RNA"
Idents(data.combined) <- "RNA_snn_res.0.9"

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/129_FindAllMarkers_HarmonySample_res0.9_MG_subcluster.tsv")

# pick the top 100 markers per cluster
sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  write_tsv("../../out/table/129_FindAllMarkers_HarmonySample_res0.9_MG_subcluster_top100.tsv")

sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  write_tsv("../../out/table/129_FindAllMarkers_HarmonySample_res0.9_MG_subcluster_top100_noRIBOandMT.tsv")

# try plotting the top markers
top_specific_markers <- sobj_total_h.markers %>%
  # filter ribosomal and mt genes
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)

# And generate e.g. a dotplot:
dittoSeq::dittoDotPlot(data.combined,
                       vars = unique(top_specific_markers$gene), 
                       group.by = "RNA_snn_res.0.9")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
ggsave("../../out/image/129_Ditto_HarmonySample_res0.9_MG_subcluster.pdf",width = 15,height = 5)

# plot the proportions
df_summary <- data.combined@meta.data %>%
  group_by(dataset,pathology_class,origin,RNA_snn_res.0.9) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(dataset,pathology_class,origin) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

df_summary %>%
  ggplot(aes(x=pathology_class,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2),shape=1)+
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1,angle = 45)) +
  facet_wrap(~RNA_snn_res.0.9,ncol = 1,scales = "free")+
  scale_y_sqrt()
ggsave("../../out/image/129_plot_clusterProp_res0.9_MG_subcluster.pdf",height = 30,width = 4)

# try the same with propeller on the same data
# renv::install("phipsonlab/speckle")
# renv::install("statmod")
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)


# run the test ------------------------------------------------------------
meta_test <- data.combined@meta.data
# filter(treat != "CSF.MS_RAPA")

table(meta_test$dataset,meta_test$orig.ident)
table(meta_test$pathology_class,meta_test$dataset)
out_diagnosis <- propeller(clusters = meta_test$RNA_snn_res.0.9,
                           sample = paste0(meta_test$dataset),
                           group = meta_test$pathology_class)

out_diagnosis %>%
  rownames_to_column("RNA_snn_res.0.9") %>%
  write_tsv("../../out/table/129_propeller_res0.9_MG_subcluster.tsv")

# plotting diagnosis ------------------------------------------------------
# default plot
speckle::plotCellTypeProps(x = data.combined,
                           clusters = data.combined$RNA_snn_res.0.9,
                           sample = data.combined$pathology_class)+theme_minimal()+theme(panel.grid = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/129_plot_propeller_res0.9_MG_subcluster.pdf",height = 8,width = 5)

# custom plot
df_summary_diagnosis <- meta_test %>% 
  mutate(group_id = dataset) %>%
  group_by(RNA_snn_res.0.9,
           group_id,
           pathology_class) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(group_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

df_summary_diagnosis %>%
  write_tsv("../../out/table/129_df_summary_diagnosis_res0.9_MG_subcluster.tsv")

# plot 01
df_summary_diagnosis %>%
  ggplot(aes(x=pathology_class,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~RNA_snn_res.0.9,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/129_propeller_res0.9_plot01_MG_subcluster.pdf",width = 10,height = 10)

# plot 02
df_summary_diagnosis %>%
  ggplot() +
  geom_boxplot(aes(x=RNA_snn_res.0.9,y=prop,color=pathology_class),outlier.shape = NA) +
  geom_point(aes(x=RNA_snn_res.0.9,y=prop,color=pathology_class),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
# ggsave("../../out/image/manualClean/propeller_plot02_diagnosis_cellid.pdf",width = 8,height = 5)

# plot expression ---------------------------------------------------------
# Martina asked to plot CD83 expression in this dataset
FeaturePlot(data.combined,features = "CD83")+scale_color_viridis_c(option = "turbo")

# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
# data.combined$group <- paste0(data.combined$pathology_class,"-",data.combined$RNA_snn_res.0.9,"-",data.combined$orig.ident)
data.combined$group <- paste0(data.combined$RNA_snn_res.0.9)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,group.by = c("group"))

GOI <- c("CD83","VEGFA","NAMPT","HIF1A","SULF2","FOSL2","REL","HK2","TNFRSF1B","DDIT4")

df_avg <- average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  filter(gene %in% GOI) %>%
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  # mutate(pathology_class = str_extract(group,pattern = c("CX_Ctrl|CX_Demye|CX_Mye|WM_CA|WM_CI|WM_Core|WM_Ctrl|WM_NAWM"))) |> 
  # mutate(donor = str_extract(group,pattern = c("s\\d+"))) |> 
  # mutate(RNA_snn_res.0.9 = str_remove_all(group,pattern = c("CX_Ctrl|CX_Demye|CX_Mye|WM_CA|WM_CI|WM_Core|WM_Ctrl|WM_NAWM|s\\d+"))%>% str_remove_all(pattern = "\\."))
  mutate(RNA_snn_res.0.9 = str_remove_all(group,pattern = c("X")))

# plot the heatmap
mat <- df_avg %>%
  group_by(gene) %>%
  mutate(scaled_avg_exp = (avg_exp - mean(avg_exp))/sd(avg_exp)) %>%
  ungroup() %>%
  select(gene,RNA_snn_res.0.9,scaled_avg_exp) %>%
  pivot_wider(names_from = RNA_snn_res.0.9,values_from = scaled_avg_exp) %>%
  column_to_rownames("gene")

Heatmap(mat,name = "scaled_avg_exp",col = colorRamp2(c(-1, 0, 2), c("blue", "white", "red")))
Heatmap(mat,name = "scaled_avg_exp",col = viridis::viridis(option = "turbo",n = 10))

# pdf("../../out/image/129_heatmap_senescence_res0.9_IMMLYMsubset_highres.pdf",width = 12,height = 5)
draw(hm03,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
# dev.off()

# # plot the average expresison by cell annotation
# df_avg |>
#   # ggplot(aes(x=NMDA_time,y=count)) + 
#   ggplot(aes(x=RNA_snn_res.0.9,y=avg_exp))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
#   # geom_col()+
#   # facet_wrap(~cell_type2,scales = "free")+
#   theme_bw()+
#   theme(axis.text.x = element_text(hjust = 1,angle = 90))+
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))+
#   facet_wrap(~gene,scales = "free")+scale_y_sqrt()
# 
# # do the same as above but split by condition
# df_avg |>
#   # ggplot(aes(x=NMDA_time,y=count)) + 
#   ggplot(aes(x=pathology_class,y=avg_exp))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
#   # geom_col()+
#   # facet_wrap(~cell_type2,scales = "free")+
#   theme_bw()+
#   theme(axis.text.x = element_text(hjust = 1,angle = 90))+
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))+
#   facet_wrap(~gene,scales = "free")
# 
# # plot splitting by treat full
# df_avg |>
#   # ggplot(aes(x=NMDA_time,y=count)) + 
#   ggplot(aes(x=pathology_class,y=avg_exp))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
#   # geom_col()+
#   # facet_wrap(~cell_type2,scales = "free")+
#   theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))+
#   facet_wrap(~RNA_snn_res.0.9,scales = "free")
# # scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# # ggsave("../../out/image/06_dotplot_annotationConfident_GLP1R_expressionAvg_treatFull.pdf",width = 9,height = 9)
# 
# # plot splitting by treat
# df_avg |>
#   # ggplot(aes(x=NMDA_time,y=count)) + 
#   ggplot(aes(x=pathology_class,y=avg_exp))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
#   # geom_col()+
#   # facet_wrap(~cell_type2,scales = "free")+
#   theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))+
#   facet_wrap(~RNA_snn_res.0.9,scales = "free")
# # scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# # ggsave("../../out/image/00_dotplot_annotationConfident_GLP1R_expressionAvg_treat.pdf",width = 9,height = 9)
# 
# # try to keep the same scale
# # calculate the median per annotation
# df_avg_summary <- df_avg %>%
#   group_by(expertAnno.l1) %>%
#   summarise(med = median(avg_exp)) %>%
#   ungroup() %>%
#   mutate(expertAnno.l1 = fct_reorder(expertAnno.l1,desc(med)))
# 
# df_avg %>%
#   mutate(expertAnno.l1 = factor(expertAnno.l1,levels = levels(df_avg_summary$expertAnno.l1))) %>%
#   # ggplot(aes(x=NMDA_time,y=count)) + 
#   ggplot(aes(x=pathology_class,y=avg_exp))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
#   geom_hline(data = df_avg_summary,aes(yintercept = med),col="red",linetype="dashed") +
#   # geom_col()+
#   # facet_wrap(~cell_type2,scales = "free")+
#   theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))+
#   facet_wrap(~expertAnno.l1,nrow=1)
# # scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/00_dotplot_annotationConfident_GLP1R_expressionAvg_treat_scale.pdf",width = 14,height = 3)