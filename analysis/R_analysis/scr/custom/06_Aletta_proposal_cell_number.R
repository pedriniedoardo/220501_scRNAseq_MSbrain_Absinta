# AIM ---------------------------------------------------------------------
# calculate the number of cells per sample and condition

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
p_umap_global <- DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")
ggsave(plot = p_umap_global,filename = "../../out/image/custom/06_UMAP_global_annotation.pdf",width = 8,height = 6)

# load the finer annation for the immune subcluster
test <- readRDS("../../out/object/129_MG_subcluster_HarmonySample.rds")
# DimPlot(test,label = T)

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
p_umap_immune <- DimPlot(test,label = T)
ggsave(plot = p_umap_immune,filename = "../../out/image/custom/06_UMAP_immune_subcluster_annotation.pdf",width = 8,height = 6)

# Idents(test) <- "cluster_martina_numeric"
# DimPlot(test,label = T)

# wrangling ---------------------------------------------------------------
# extract the meta from the full object
df_01 <- data.combined@meta.data %>%
  group_by(orig.ident) %>%
  summarise(n_tot = n())


df_02 <- data.combined@meta.data %>%
  group_by(orig.ident,sample_id,facility,origin,disease,sex,age,pathology,pathology_class,patient,plaque,origin_alt,expertAnno.l1) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = expertAnno.l1,
              values_from = n,
              names_glue = "n_{expertAnno.l1}",
              values_fill = 0)

df_full <- left_join(df_02,df_01,by = "orig.ident")

# extract the meta from the MG subcluster using the re-annotation suggested by Martina
df_IMM_01 <- test@meta.data %>%
  group_by(orig.ident) %>%
  summarise(n_tot = n())


df_IMM_02 <- test@meta.data %>%
  group_by(orig.ident,sample_id,facility,origin,disease,sex,age,pathology,pathology_class,patient,plaque,origin_alt,cluster_martina) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = cluster_martina,
              values_from = n,
              names_glue = "n_{cluster_martina}",
              values_fill = 0)

df_IMM_full <- left_join(df_IMM_02,df_IMM_01,by = "orig.ident")

# saving --------------------------------------------------------------------
write_tsv(df_full,"../../out/table/custom/06_Aletta_cellNumber_global.tsv")
write_tsv(df_IMM_full,"../../out/table/custom/06_Aletta_cellNumber_immune_subcluster.tsv")

