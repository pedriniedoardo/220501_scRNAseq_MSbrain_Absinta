# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(dittoSeq)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v3")

# -------------------------------------------------------------------------
test <- readRDS("../../out/object/129_MG_subcluster_HarmonySample.rds")

DimPlot(test,label = T)

# martina suggested to regreoup the cells following the current clustering
meta_new <- test@meta.data %>%
  mutate(cluster_martina = case_when(RNA_snn_res.0.9 %in% c(18) ~ "plasmablast",
                                     RNA_snn_res.0.9 %in% c(8) ~ "T-cells",
                                     RNA_snn_res.0.9 %in% c(10) ~ "per-mac",
                                     RNA_snn_res.0.9 %in% c(6,4) ~ "dendritic",
                                     RNA_snn_res.0.9 %in% c(12,19,5) ~ "mye-MG",
                                     RNA_snn_res.0.9 %in% c(1,11) ~ "cont",
                                     RNA_snn_res.0.9 %in% c(7,16) ~ "MIMS-foamy",
                                     RNA_snn_res.0.9 %in% c(9) ~ "MIMS-iron",
                                     RNA_snn_res.0.9 %in% c(13) ~ "stress-MG",
                                     RNA_snn_res.0.9 %in% c(0,3) ~ "homeo-MG",
                                     RNA_snn_res.0.9 %in% c(15) ~ "cl9Absinta",
                                     RNA_snn_res.0.9 %in% c(17) ~ "cl17",
                                     RNA_snn_res.0.9 %in% c(14) ~ "cl14",
                                     RNA_snn_res.0.9 %in% c(2) ~ "und"))

# define the order of the factors
order_factors <- meta_new %>%
  group_by(cluster_martina) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  pull(cluster_martina)

meta_new_final <- meta_new %>%
  mutate(cluster_martina_numeric = factor(cluster_martina,levels = order_factors) %>% as.numeric()) %>%
  mutate(cluster_martina_numeric = cluster_martina_numeric - 1) %>%
  mutate(cluster_martina_numeric = as.factor(cluster_martina_numeric))

# add the metadata to the original object
test$cluster_martina_numeric <- meta_new_final$cluster_martina_numeric
test$cluster_martina <- meta_new_final$cluster_martina

Idents(test) <- "cluster_martina"
DimPlot(test,label = T)
ggsave("../../out/image/01_MG_subset_clustersMartina.pdf",width = 7,height = 6)

Idents(test) <- "cluster_martina_numeric"
DimPlot(test,label = T)
ggsave("../../out/image/01_MG_subset_clustersMartinaNumeric.pdf",width = 5,height = 4)

# save the object with the new annotation
saveRDS(test,"../../out/object/129_MG_subcluster_HarmonySample_martinaCluster.rds")

# make the marker genes per cluster
# Idents(test)
test.markers <- FindAllMarkers(test, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_tsv(test.markers,"../../out/table/01_MG_subset_clustersMartina_markers.tsv")

# pick the top 100 markers per cluster
test.markers %>%
  group_by(cluster) %>%
  slice(1:100) %>%
  write_tsv("../../out/table/01_MG_subset_clustersMartina_markersTop100.tsv")

# remove the technical genes
test.markers %>%
  group_by(cluster) %>%
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  slice(1:100) %>%
  write_tsv("../../out/table/01_MG_subset_clustersMartina_markersTop100_noRIBOandMT.tsv")

# try plotting the top markers
top_markers <- test.markers %>%
  # filter ribosomal and mt genes
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC)

# And generate e.g. a dotplot:
dittoSeq::dittoDotPlot(test,
                       vars = unique(top_markers$gene), 
                       group.by = "cluster_martina_numeric")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
ggsave("../../out/image/01_MG_subset_clustersMartina_dittoseq.pdf",width = 12,height = 4)

# make the dotplot using the panel
# martina asked to remove FTL
shortlist_features_list <- list(
  IMMUNE = c("LYVE1", "CD163", "MRC1", "LINGO1", "HSPA1A", "MOBP", "CD83", "HIF1A", "CD22", "VEGFA", "SOD1", "TREM2", "CX3CR1", "P2RY12", "C3", "CSF1R", "CD74", "RUNX1", "C1QB", "PTPRC", "AIF1", "HLA-DRA", "TYROBP"),
  B_CELLS = c("IGHG1", "CD38"),
  T_CELLS = c("SKAP1", "CD8A", "CD2")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
Idents(test) <- "cluster_martina_numeric"
dotplot <- DotPlot(test, features = shortlist_features_list, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
# ggsave(plot = dotplot,"../out/plot/01_MG_subset_clustersMartina_dotplot.pdf",width = 13,height = 5)

# try ggplot implementation
df_test <- lapply(shortlist_features_list,function(x){
  dotplot$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

test_long01 <- df_test %>%
  # force the order
  # mutate(id = factor(id,levels = c(9,5,6,3,2,7,4,8,0,1,10))) %>% 
  mutate(cell_type = factor(cell_type,levels = c("IMMUNE","B_CELLS","T_CELLS"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 45),
        strip.text.x = element_text(angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue",limits = c(-1,2),oob = scales::squish)
ggsave(plot = test_long01,"../../out/image/01_MG_subset_clustersMartina_dotplot.pdf",width = 13,height = 5)

# -------------------------------------------------------------------------
# add module score assessment
df_LUT <- data.frame(type = c("homeostatic","MIMs foamy","MIMs iron","DC_2","DC_6","stressed MG"),
                     cluster = c(0,1,8,2,6,3))

list_sig <- read_csv("../../data/Signature_IMM_subset_Absinta2019.csv") %>%
  left_join(df_LUT,by = "cluster") %>%
  mutate(type = case_when(is.na(type)~paste0("clu_",cluster),
                          T~type)) %>%
  split(.$type) %>%
  lapply(function(x){
    x %>%
      pull(gene) %>%
      unique()
  })

# score the signatures for senescence
test <- Seurat::AddModuleScore(test,
                               features = list_sig,
                               name = "_score")

df_rename_long <- data.frame(names = test@meta.data %>%
                               colnames() %>%
                               str_subset("_score"),
                             rename = paste0("score_",names(list_sig)))

lookup_long <- df_rename_long$names
names(lookup_long) <- df_rename_long$rename

# rename the columns
test@meta.data <- dplyr::rename(test@meta.data,all_of(lookup_long))

# plot the scores from AddModuleScore
list_plot_02_long <- lapply(df_rename_long$rename,function(x){
  plot <- FeaturePlot(test,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_long)
ggsave("../../out/image/02_UMAP_SigAbsinta2019_IMMLYMsubset_highres.pdf",width = 25,height = 20)

# plot the average score per signature per cluster as a heatmap
df_senescence <- test@meta.data %>%
  dplyr::select(cluster_martina_numeric,contains("score_")) %>%
  pivot_longer(names_to = "signature",values_to = "score",-cluster_martina_numeric) %>%
  group_by(cluster_martina_numeric,signature) %>%
  summarise(avg_score = mean(score),
            med_score = median(score)) %>%
  mutate(cluster_id = paste0("clu_",cluster_martina_numeric,"_custom")) %>%
  ungroup()

mat_senescence_avg <- df_senescence %>%
  group_by(signature) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = cluster_id,values_from = scaled_avg_score) %>%
  column_to_rownames("signature")

hm01 <- Heatmap(mat_senescence_avg,name = "scaled_avg_exp",col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")))

mat_senescence_avg %>%
  rownames_to_column("signatures") %>%
  write_tsv("../../out/table/02_tableScaled_SigAbisnta2019_custom_IMMLYMsubset_highres.tsv")

df_senescence %>%
  write_tsv("../../out/table/02_tableUnScaled_SigAbisnta2019_custom_IMMLYMsubset_highres.tsv")

pdf("../../out/image/02_heatmap_SigAbisnta2019_custom_IMMLYMsubset_highres.pdf",width = 9,height = 4)
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

pdf("../../out/image/02_heatmap_SigAbisnta2019_custom_IMMLYMsubset_highres2.pdf",width = 8,height = 5)
draw(hm02,heatmap_legend_side = "left",padding = unit(c(4, 2, 2, 80), "mm"))
dev.off()
