# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration. in particular to plot the senescence cells

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
data.combined <- readRDS("../../out/object/129_MG_subcluster_HarmonySample.rds")
# data.combined2 <- readRDS("../../out/object/100_MG_subcluster_HarmonyRun.rds")
(DimPlot(data.combined,group.by = "dataset") + ggtitle("Harmony Sample"))
(DimPlot(data.combined,group.by = "pathology_class",split.by = "dataset") + ggtitle("Harmony Sample"))
(DimPlot(data.combined,group.by = "RNA_snn_res.0.9") + ggtitle("Harmony Sample"))

# pull the annotation derived from the SIT algorithm
df_SIT <- read_tsv("../../out/table/revision/124_meta_SIT_WM_CX_harmonySkipIntegration_AllSoupX_test.tsv")

# pull the annotation derived from SENMAYO
df_SENMAYO <- read_tsv("../../out/table/121_meta_SENMAYO_WM_CX_harmonySkipIntegration_AllSoupX_test.tsv")

# wrangling ---------------------------------------------------------------
# assign the call of senescence to each cell in the subcluster analysis
df_meta_full <- data.combined@meta.data %>%
  rownames_to_column("barcodes") %>%
  left_join(df_SIT %>% dplyr::select(barcodes = rowname,sen_SIT = SENEQUANTILE) %>% mutate(sen_SIT = case_when(sen_SIT == "YES"~1,T~0))) %>%
  left_join(df_SENMAYO %>% dplyr::select(barcodes = barcode,sen_SENMAYO = sen)) %>%
  column_to_rownames("barcodes")

table(df_meta_full$sen_SIT,df_meta_full$sen_SENMAYO)

# swap the metadata
data.combined@meta.data <- df_meta_full


# define the proportion per cluster using the two methods
# regardless of the sample
df_meta_full %>%
  group_by(RNA_snn_res.0.9) %>%
  summarise(n = n(),
            tot_sen_SIT = sum(sen_SIT),
            tot_sen_SENMAYO = sum(sen_SENMAYO)) %>%
  ungroup() %>%
  mutate(prop_SIT = tot_sen_SIT/n,
         prop_SENMAYO = tot_sen_SENMAYO/n) %>%
  pivot_longer(names_to = "var",values_to = "prop",c(prop_SIT,prop_SENMAYO)) %>%
  group_by(RNA_snn_res.0.9) %>%
  mutate(avg_prop = mean(prop)) %>%
  ungroup() %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,-avg_prop)) %>%
  ggplot(aes(x=RNA_snn_res.0.9,y=prop,fill=var))+geom_col(position = "dodge")+theme_cowplot()

# regardless of the sample
df_meta_full %>%
  group_by(RNA_snn_res.0.9) %>%
  summarise(n = n(),
            tot_sen_SIT = sum(sen_SIT),
            tot_sen_SENMAYO = sum(sen_SENMAYO)) %>%
  ungroup() %>%
  mutate(prop_SIT = tot_sen_SIT/n,
         prop_SENMAYO = tot_sen_SENMAYO/n) %>%
  pivot_longer(names_to = "var",values_to = "prop",c(prop_SIT,prop_SENMAYO)) %>%
  group_by(RNA_snn_res.0.9) %>%
  mutate(avg_prop = mean(prop)) %>%
  ungroup() %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,-avg_prop)) %>%
  ggplot(aes(x=RNA_snn_res.0.9,y=prop,fill=var))+geom_col(position = "dodge")+theme_cowplot()

# save the ref table
df_meta_full %>%
  group_by(RNA_snn_res.0.9) %>%
  summarise(n = n(),
            tot_sen_SIT = sum(sen_SIT),
            tot_sen_SENMAYO = sum(sen_SENMAYO)) %>%
  ungroup() %>%
  mutate(prop_SIT = tot_sen_SIT/n,
         prop_SENMAYO = tot_sen_SENMAYO/n) %>%
  pivot_longer(names_to = "var",values_to = "prop",c(prop_SIT,prop_SENMAYO)) %>%
  group_by(RNA_snn_res.0.9) %>%
  mutate(avg_prop = mean(prop)) %>%
  ungroup() %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,-avg_prop)) %>%
  write_tsv("../../out/table/129_prop_res0.9_senescence.tsv")

# c(1,11,14,5,12,19,8,15,17,18)
df_meta_full %>%
  filter(!RNA_snn_res.0.9 %in% c(1,11,14,5,12,19,8,15,17,18)) %>%
  group_by(RNA_snn_res.0.9) %>%
  summarise(n = n(),
            tot_sen_SIT = sum(sen_SIT),
            tot_sen_SENMAYO = sum(sen_SENMAYO)) %>%
  ungroup() %>%
  mutate(prop_SIT = tot_sen_SIT/n,
         prop_SENMAYO = tot_sen_SENMAYO/n) %>%
  pivot_longer(names_to = "var",values_to = "prop",c(prop_SIT,prop_SENMAYO)) %>%
  group_by(RNA_snn_res.0.9) %>%
  mutate(avg_prop = mean(prop)) %>%
  ungroup() %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,-avg_prop)) %>%
  ggplot(aes(x=RNA_snn_res.0.9,y=prop,fill=var))+geom_col(position = "dodge")+theme_cowplot()

# add the sample info
df_meta_full %>%
  filter(!RNA_snn_res.0.9 %in% c(1,11,14,5,12,19,8,15,17,18)) %>%
  group_by(orig.ident,pathology_class,RNA_snn_res.0.9) %>%
  summarise(n = n(),
            tot_sen_SIT = sum(sen_SIT),
            tot_sen_SENMAYO = sum(sen_SENMAYO)) %>%
  ungroup() %>%
  mutate(prop_SIT = tot_sen_SIT/n,
         prop_SENMAYO = tot_sen_SENMAYO/n) %>%
  pivot_longer(names_to = "var",values_to = "prop",c(prop_SIT,prop_SENMAYO)) %>%
  group_by(RNA_snn_res.0.9) %>%
  mutate(avg_prop = mean(prop)) %>%
  ungroup() %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,-avg_prop)) %>%
  ggplot(aes(x=RNA_snn_res.0.9,y=prop)) +
  geom_boxplot(position = "dodge",outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),alpha=0.6) +
  theme_cowplot()+
  facet_wrap(~var,nrow=2)+
  theme(strip.background = element_blank())

# split by disease
df_meta_full %>%
  filter(!RNA_snn_res.0.9 %in% c(1,11,14,5,12,19,8,15,17,18)) %>%
  group_by(orig.ident,pathology_class,RNA_snn_res.0.9) %>%
  summarise(n = n(),
            tot_sen_SIT = sum(sen_SIT),
            tot_sen_SENMAYO = sum(sen_SENMAYO)) %>%
  ungroup() %>%
  mutate(prop_SIT = tot_sen_SIT/n,
         prop_SENMAYO = tot_sen_SENMAYO/n) %>%
  pivot_longer(names_to = "var",values_to = "prop",c(prop_SIT,prop_SENMAYO)) %>%
  group_by(RNA_snn_res.0.9) %>%
  mutate(avg_prop = mean(prop)) %>%
  ungroup() %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,-avg_prop)) %>%
  ggplot(aes(x=RNA_snn_res.0.9,y=prop,col=var)) +
  geom_boxplot(position = "dodge",outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),alpha=0.6) +
  theme_cowplot()+
  facet_wrap(~pathology_class,nrow=2)+
  theme(strip.background = element_blank())

# save the table as ref
df_meta_full %>%
  filter(!RNA_snn_res.0.9 %in% c(1,11,14,5,12,19,8,15,17,18)) %>%
  group_by(orig.ident,pathology_class,RNA_snn_res.0.9) %>%
  summarise(n = n(),
            tot_sen_SIT = sum(sen_SIT),
            tot_sen_SENMAYO = sum(sen_SENMAYO)) %>%
  ungroup() %>%
  mutate(prop_SIT = tot_sen_SIT/n,
         prop_SENMAYO = tot_sen_SENMAYO/n) %>%
  pivot_longer(names_to = "var",values_to = "prop",c(prop_SIT,prop_SENMAYO)) %>%
  group_by(RNA_snn_res.0.9) %>%
  mutate(avg_prop = mean(prop)) %>%
  ungroup() %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,-avg_prop)) %>%
  write_tsv("../../out/table/129_prop_res0.9_senescence_sample.tsv")


# check the overall trend of the cells in each cluster
df_meta_full %>%
  group_by(orig.ident,pathology_class,RNA_snn_res.0.9) %>%
  summarise(n_cluster = n()) %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(n_tot_sample = sum(n_cluster)) %>%
  ungroup() %>%
  mutate(prop = n_cluster/n_tot_sample) %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,prop)) %>%
  ggplot(aes(x=pathology_class,y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),alpha=0.6) +
  theme_cowplot()+
  facet_wrap(~RNA_snn_res.0.9,scales = "free")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))

# save the ref table
df_meta_full %>%
  group_by(orig.ident,pathology_class,RNA_snn_res.0.9) %>%
  summarise(n_cluster = n()) %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(n_tot_sample = sum(n_cluster)) %>%
  ungroup() %>%
  mutate(prop = n_cluster/n_tot_sample) %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,prop)) %>%
  write_tsv("../../out/table/129_prop_res0.9.tsv")

# split by origin and disease per cluster
df_meta_full %>%
  group_by(orig.ident,origin,disease,RNA_snn_res.0.9) %>%
  summarise(n_cluster = n()) %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(n_tot_sample = sum(n_cluster)) %>%
  ungroup() %>%
  mutate(prop = n_cluster/n_tot_sample) %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,prop)) %>%
  ggplot(aes(x=disease,y=prop,color=origin)) +
  geom_boxplot(outlier.shape = NA,position = "dodge") +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),alpha=0.6) +
  theme_cowplot()+
  facet_wrap(~RNA_snn_res.0.9,scales = "free")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))

# save the ref table
df_meta_full %>%
  group_by(orig.ident,origin,disease,RNA_snn_res.0.9) %>%
  summarise(n_cluster = n()) %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(n_tot_sample = sum(n_cluster)) %>%
  ungroup() %>%
  mutate(prop = n_cluster/n_tot_sample) %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,prop)) %>%
  write_tsv("../../out/table/129_prop_res0.9_originDisease.tsv")

# martina also wanted to include the zeros for the missing combinations
LUT_sample <- df_meta_full %>%
  group_by(orig.ident,origin,disease,pathology_class) %>%
  summarise(.groups = "drop")

df_meta_full %>%
  group_by(orig.ident,pathology_class,RNA_snn_res.0.9) %>%
  summarise(n_cluster = n(),.groups = "drop") %>%
  complete(orig.ident,RNA_snn_res.0.9, fill = list(n_cluster = 0)) %>%
  select(-pathology_class) %>%
  left_join(LUT_sample,by = "orig.ident") %>%
  group_by(orig.ident) %>%
  mutate(n_tot_sample = sum(n_cluster)) %>%
  ungroup() %>%
  mutate(prop = n_cluster/n_tot_sample) %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,prop)) %>%
  ggplot(aes(x=pathology_class,y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),alpha=0.6) +
  theme_cowplot()+
  facet_wrap(~RNA_snn_res.0.9,scales = "free")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))

df_meta_full %>%
  group_by(orig.ident,pathology_class,RNA_snn_res.0.9) %>%
  summarise(n_cluster = n(),.groups = "drop") %>%
  complete(orig.ident,RNA_snn_res.0.9, fill = list(n_cluster = 0)) %>%
  select(-pathology_class) %>%
  left_join(LUT_sample,by = "orig.ident") %>%
  group_by(orig.ident) %>%
  mutate(n_tot_sample = sum(n_cluster)) %>%
  ungroup() %>%
  mutate(prop = n_cluster/n_tot_sample) %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,prop)) %>%
  write_tsv("../../out/table/129_prop_res0.9_includeZeros.tsv")

# split by origin and disease per cluster
df_meta_full %>%
  group_by(orig.ident,origin,disease,RNA_snn_res.0.9) %>%
  summarise(n_cluster = n(),.groups = "drop") %>%
  complete(orig.ident,RNA_snn_res.0.9, fill = list(n_cluster = 0)) %>%
  select(-c(origin,disease)) %>%
  left_join(LUT_sample,by = "orig.ident") %>%
  group_by(orig.ident) %>%
  mutate(n_tot_sample = sum(n_cluster)) %>%
  ungroup() %>%
  mutate(prop = n_cluster/n_tot_sample) %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,prop)) %>%
  ggplot(aes(x=disease,y=prop,color=origin)) +
  geom_boxplot(outlier.shape = NA,position = "dodge") +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),alpha=0.6) +
  theme_cowplot()+
  facet_wrap(~RNA_snn_res.0.9,scales = "free")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))

# save the ref table
df_meta_full %>%
  group_by(orig.ident,origin,disease,RNA_snn_res.0.9) %>%
  summarise(n_cluster = n(),.groups = "drop") %>%
  complete(orig.ident,RNA_snn_res.0.9, fill = list(n_cluster = 0)) %>%
  select(-c(origin,disease)) %>%
  left_join(LUT_sample,by = "orig.ident") %>%
  group_by(orig.ident) %>%
  mutate(n_tot_sample = sum(n_cluster)) %>%
  ungroup() %>%
  mutate(prop = n_cluster/n_tot_sample) %>%
  mutate(RNA_snn_res.0.9 = fct_reorder(RNA_snn_res.0.9,prop)) %>%
  write_tsv("../../out/table/129_prop_res0.9_originDisease_includeZeros.tsv")
