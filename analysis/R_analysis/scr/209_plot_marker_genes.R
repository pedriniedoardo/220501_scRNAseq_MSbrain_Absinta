# AIM ---------------------------------------------------------------------
# aggregate and merge data to make a dotplot of know marker genes

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(DESeq2)
library(cowplot)

# read in the data --------------------------------------------------------
# read in the sample seurat object
scobj <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
DimPlot(scobj,raster = T,group.by = "expertAnno.l1",label = T)

# define a panel fo marker genes
marker_gene_panel <- list(
  IMMUNE = c("AIF1", "CX3CR1", "P2RY12", "TMEM119", "CSF1R"),
  B_CELLS = c("MS4A1", "CD79A", "JCHAIN"),
  T_CELLS = c("CD3D", "CD3E", "CD2", "CD8A"),
  OLIGO = c("MOG", "MBP", "PLP1", "MAG", "CNP"),
  OPC = c("PDGFRA", "SOX10", "OLIG2", "OLIG1"),
  ASTRO = c("GFAP", "AQP4", "ALDH1L1", "S100B", "SLC1A2"),
  NEURONS = c("RBFOX3", "SYP", "SNAP25", "SYT1", "TUBB3"),
  ENDOTHELIAL = c("CLDN5", "PECAM1", "VWF", "CDH5", "FLT1"),
  PERICYTE = c("PDGFRB", "RGS5", "ACTA2")
)

# make the list as a data.frame()
shortlist_features_df <- lapply(marker_gene_panel,function(x){
  data.frame(gene = x)
}) %>%
  bind_rows(.id = "cat")

# get the average expression per sample
# set the idents
Idents(scobj) <- "expertAnno.l1"
# DefaultAssay(data.combined) <- "RNA"

# scale the gene expression within genes. highlight the samples were each gene is mostly expressed
df_average <- AverageExpression(scobj,features = unique(shortlist_features_df$gene),slot = "data")$RNA %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene) %>%
  group_by(gene) %>%
  mutate(avg_exp.scaled = scale(avg_exp)[,1]) %>%
  # add the gene category
  left_join(shortlist_features_df,by = "gene") %>%
  ungroup()

# confirm the scaling
df_average %>%
  group_by(gene) %>%
  summarise(avg = mean(avg_exp.scaled),
            sd = sd(avg_exp.scaled)) %>%
  print(n=40)

# make ggplot heatmap
df_average %>%
  # force the order
  mutate(group_id = factor(group_id,levels = c("AST","LYM","VAS","IMM","EXC NEU","INH NEU","OLIGO","OPC","EPENDYMA"))) %>% 
  mutate(cat = factor(cat,levels = c("ASTRO","B_CELLS","T_CELLS","ENDOTHELIAL","PERICYTE","IMMUNE","NEURONS","OLIGO","OPC"))) %>% 
  ggplot(aes(x = gene,y = group_id)) +
  # geom_tile(aes(fill = avg_exp.scaled))+
  geom_tile(aes(fill = avg_exp.scaled), width = 0.95, height = 0.95, size = 0.2, color = "white")+
  # scale_radius(range = c(0, 8)) +
  facet_grid(~cat,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 0))+
  scale_fill_gradientn(colours = viridis::turbo(20), limits = c(-0.5,3), oob = scales::squish, name = 'avg exp scaled\n within genes')
ggsave("../../out/image/209_heatmap_scale_within_genes_sc.pdf",width = 14,height = 4)

# scale the gene expression between genes. highlight the most expressed genes
df_average2 <- AverageExpression(scobj,features = unique(shortlist_features_df$gene),slot = "data")$RNA %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene) %>%
  group_by(group_id) %>%
  mutate(avg_exp.scaled = scale(avg_exp)[,1]) %>%
  # add the gene category
  left_join(shortlist_features_df,by = "gene") %>%
  ungroup()

# confirm the scaling
df_average2 %>%
  group_by(group_id) %>%
  summarise(avg = mean(avg_exp.scaled),
            sd = sd(avg_exp.scaled)) %>%
  print(n=40)

# make ggplot heatmap
df_average2 %>%
  # force the order
  mutate(group_id = factor(group_id,levels = c("AST","LYM","VAS","IMM","EXC NEU","INH NEU","OLIGO","OPC","EPENDYMA"))) %>% 
  mutate(cat = factor(cat,levels = c("ASTRO","B_CELLS","T_CELLS","ENDOTHELIAL","PERICYTE","IMMUNE","NEURONS","OLIGO","OPC"))) %>% 
  ggplot(aes(x = gene,y = group_id)) +
  # geom_tile(aes(fill = avg_exp.scaled))+
  geom_tile(aes(fill = avg_exp.scaled), width = 0.95, height = 0.95, size = 0.2, color = "white")+
  # scale_radius(range = c(0, 8)) +
  facet_grid(~cat,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 0))+
  scale_fill_gradientn(colours = viridis::turbo(20),
                       # limits = c(0,2),
                       oob = scales::squish, name = 'avg exp scaled\n between genes')
ggsave("../../out/image/209_heatmap_scale_between_genes_sc.pdf",width = 14,height = 4)

# do the same on the bulk data from Eliana --------------------------------
# read Eliana's bulk data red the input before filtering low expressed genes
# dds <- readRDS("../../out/object/201_dds_all.rds") %>%
#   DESeq()
dds <- readRDS("../../out/object/201_dds_all_update.rds") %>%
  DESeq()

# extract the normalized table of counts
exp_norm <- counts(dds,normalized = T) %>%
  data.frame() %>%
  # dplyr::select(contains("BASELINE")|contains("Fe")|contains("myelin")) %>%
  as.matrix()

# make it a long format
exp_norm_long <- exp_norm %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "sample_id",values_to = "norm_exp",-gene)

# join the tables, if a gene is missing coerce the expression to 0
# build a reference matrix to be filled. I should expect each gene is present in each sample
# This assumes that a gene not being detected is equivalent to it having zero expression. While practical, this is a strong assumption. A gene might be absent due to low expression falling below the detection threshold, not true absence

df_average_bulk <- crossing(gene = unique(shortlist_features_df$gene),
         sample_id = colnames(exp_norm)) %>%
  left_join(exp_norm_long,by = c("gene","sample_id")) %>%
  # filter(is.na(norm_exp))
  mutate(norm_exp = case_when(is.na(norm_exp)~0,
                                    T~norm_exp)) %>%
  # scale the expression
  group_by(gene) %>%
  mutate(norm_exp.scaled = scale(norm_exp)[,1]) %>%
  # add the gene category
  left_join(shortlist_features_df,by = "gene") %>%
  ungroup()

# confirm the scaling
df_average_bulk %>%
  group_by(gene) %>%
  summarise(avg = mean(norm_exp.scaled),
            sd = sd(norm_exp.scaled)) %>%
  print(n=40)

# make ggplot heatmap
df_average_bulk %>%
  # force the order
  # mutate(group_id = factor(group_id,levels = c("AST","LYM","VAS","IMM","EXC NEU","INH NEU","OLIGO","OPC","EPENDYMA"))) %>% 
  mutate(cat = factor(cat,levels = c("ASTRO","B_CELLS","T_CELLS","ENDOTHELIAL","PERICYTE","IMMUNE","NEURONS","OLIGO","OPC"))) %>% 
  ggplot(aes(x = gene,y = sample_id)) +
  # geom_tile(aes(fill = avg_exp.scaled))+
  geom_tile(aes(fill = norm_exp.scaled), width = 0.95, height = 0.95, size = 0.2, color = "white")+
  # scale_radius(range = c(0, 8)) +
  facet_grid(~cat,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 0))+
  scale_fill_gradientn(colours = viridis::turbo(20), limits = c(-0.5,3), oob = scales::squish, name = 'avg exp scaled\n within genes')
ggsave("../../out/image/209_heatmap_scale_within_genes_bulk.pdf",width = 14,height = 4)

# test scale across samples
df_average_bulk2 <- crossing(gene = unique(shortlist_features_df$gene),
                             sample_id = colnames(exp_norm)) %>%
  left_join(exp_norm_long,by = c("gene","sample_id")) %>%
  # filter(is.na(norm_exp))
  mutate(norm_exp = case_when(is.na(norm_exp)~0,
                              T~norm_exp)) %>%
  # scale the expression
  group_by(sample_id) %>%
  mutate(norm_exp.scaled = scale(norm_exp)[,1]) %>%
  # add the gene category
  left_join(shortlist_features_df,by = "gene") %>%
  ungroup()

# confirm the scaling
df_average_bulk2 %>%
  group_by(sample_id) %>%
  summarise(avg = mean(norm_exp.scaled),
            sd = sd(norm_exp.scaled)) %>%
  print(n=40)

df_average_bulk2 %>%
  # force the order
  # mutate(group_id = factor(group_id,levels = c("AST","LYM","VAS","IMM","EXC NEU","INH NEU","OLIGO","OPC","EPENDYMA"))) %>% 
  mutate(cat = factor(cat,levels = c("ASTRO","B_CELLS","T_CELLS","ENDOTHELIAL","PERICYTE","IMMUNE","NEURONS","OLIGO","OPC"))) %>% 
  ggplot(aes(x = gene,y = sample_id)) +
  # geom_tile(aes(fill = avg_exp.scaled))+
  geom_tile(aes(fill = norm_exp.scaled), width = 0.95, height = 0.95, size = 0.2, color = "white")+
  # scale_radius(range = c(0, 8)) +
  facet_grid(~cat,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 0))+
  scale_fill_gradientn(colours = viridis::turbo(20), limits = c(-0.5,3), oob = scales::squish, name = 'avg exp scaled\n between genes')
ggsave("../../out/image/209_heatmap_scale_between_genes_bulk.pdf",width = 14,height = 4)

# try plotting vst scaled data --------------------------------------------
vst <- vst(dds,blind = T)

# extract the normalized table of counts
exp_vst <-  assay(vst) %>%
  data.frame() %>%
  # dplyr::select(contains("BASELINE")|contains("Fe")|contains("myelin")) %>%
  as.matrix()

# make it a long format
exp_vst_long <- exp_vst %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "sample_id",values_to = "vst_exp",-gene)

# join the tables, if a gene is missing coerce the expression to 0
# build a reference matrix to be filled. I should expect each gene is present in each sample
# This assumes that a gene not being detected is equivalent to it having zero expression. While practical, this is a strong assumption. A gene might be absent due to low expression falling below the detection threshold, not true absence

df_averagevst_bulk <- crossing(gene = unique(shortlist_features_df$gene),
                            sample_id = colnames(exp_vst)) %>%
  left_join(exp_vst_long,by = c("gene","sample_id")) %>%
  # filter(is.na(norm_exp))
  mutate(norm_exp = case_when(is.na(vst_exp)~NA,
                              T~vst_exp)) %>%
  # scale the expression
  group_by(gene) %>%
  mutate(vst_exp.scaled = scale(vst_exp)[,1]) %>%
  # add the gene category
  left_join(shortlist_features_df,by = "gene") %>%
  ungroup()

# confirm the scaling
df_averagevst_bulk %>%
  group_by(gene) %>%
  summarise(avg = mean(vst_exp.scaled),
            sd = sd(vst_exp.scaled)) %>%
  print(n=40)

# make ggplot heatmap
df_averagevst_bulk %>%
  # force the order
  # mutate(group_id = factor(group_id,levels = c("AST","LYM","VAS","IMM","EXC NEU","INH NEU","OLIGO","OPC","EPENDYMA"))) %>% 
  mutate(cat = factor(cat,levels = c("ASTRO","B_CELLS","T_CELLS","ENDOTHELIAL","PERICYTE","IMMUNE","NEURONS","OLIGO","OPC"))) %>% 
  ggplot(aes(x = gene,y = sample_id)) +
  # geom_tile(aes(fill = avg_exp.scaled))+
  geom_tile(aes(fill = vst_exp.scaled), width = 0.95, height = 0.95, size = 0.2, color = "white")+
  # scale_radius(range = c(0, 8)) +
  facet_grid(~cat,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 0))+
  scale_fill_gradientn(colours = viridis::turbo(20),
                       # limits = c(-0.5,3),
                       oob = scales::squish, name = 'avg exp scaled\n within genes')
# ggsave("../../out/image/209_heatmap_scale_within_genes_bulk.pdf",width = 14,height = 4)

# test scale across samples
df_averagevst_bulk2 <- crossing(gene = unique(shortlist_features_df$gene),
                             sample_id = colnames(exp_vst)) %>%
  left_join(exp_vst_long,by = c("gene","sample_id")) %>%
  # filter(is.na(norm_exp))
  mutate(norm_exp = case_when(is.na(vst_exp)~0,
                              T~vst_exp)) %>%
  # scale the expression
  group_by(sample_id) %>%
  mutate(vst_exp.scaled = scale(vst_exp)[,1]) %>%
  # add the gene category
  left_join(shortlist_features_df,by = "gene") %>%
  ungroup()

# confirm the scaling
df_averagevst_bulk2 %>%
  group_by(sample_id) %>%
  summarise(avg = mean(vst_exp.scaled),
            sd = sd(vst_exp.scaled)) %>%
  print(n=40)

df_averagevst_bulk2 %>%
  # force the order
  # mutate(group_id = factor(group_id,levels = c("AST","LYM","VAS","IMM","EXC NEU","INH NEU","OLIGO","OPC","EPENDYMA"))) %>% 
  mutate(cat = factor(cat,levels = c("ASTRO","B_CELLS","T_CELLS","ENDOTHELIAL","PERICYTE","IMMUNE","NEURONS","OLIGO","OPC"))) %>% 
  ggplot(aes(x = gene,y = sample_id)) +
  # geom_tile(aes(fill = avg_exp.scaled))+
  geom_tile(aes(fill = vst_exp.scaled), width = 0.95, height = 0.95, size = 0.2, color = "white")+
  # scale_radius(range = c(0, 8)) +
  facet_grid(~cat,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 0))+
  scale_fill_gradientn(colours = viridis::turbo(20), limits = c(-0.5,3), oob = scales::squish, name = 'avg exp scaled\n between genes')
# ggsave("../../out/image/209_heatmap_scale_between_genes_bulk.pdf",width = 14,height = 4)



# -------------------------------------------------------------------------
# # make a regular heatmap
# mat_average <- df_average %>%
#   select(-c(avg_exp,cat)) %>%
#   pivot_wider(names_from = group_id,values_from = "avg_exp.scaled") %>%
#   column_to_rownames("gene")
# 
# Heatmap(mat_average,
#         col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
# 
# 
# 
# df_test <- lapply(shortlist_features_list_short2,function(x){
#   test_short062$data %>% 
#     filter(features.plot %in% x)
# }) %>% 
#   bind_rows(.id = "cell_type")
# 
# glimpse(df_test)
# 
# 
# # -------------------------------------------------------------------------
# 
# 
# # source image
# test_short <- DotPlot(cts_sample_CellidSmall,
#                       features = shortlist_features_list_short,
#                       dot.scale = 8,
#                       cluster.idents = T,
#                       group.by = "cell_id") +
#   RotatedAxis() +
#   labs(title = "cell_id") +
#   theme(strip.text = element_text(angle = 90))
# 
# test_short$data
# 
# 
# test_long01 <- df_test %>%
#   # force the order
#   mutate(id = factor(id,levels = c(9,5,6,3,2,7,4,8,0,1,10))) %>% 
#   mutate(cell_type = factor(cell_type,levels = c("Bc","Tc","RBc","FIBRO","MAC","DC","ENDO"))) %>% 
#   ggplot(aes(x = features.plot,y = id)) +
#   geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
#   scale_radius(range = c(0, 8)) +
#   facet_grid(~cell_type,scales = "free",space = "free")+
#   theme_cowplot()+
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(hjust = 1,angle = 90),
#         strip.text.x = element_text(angle = 0))+
#   scale_color_gradient(low = "lightgrey",high = "blue",limits = c(-1,2),oob = scales::squish)
# 
# 
# # -------------------------------------------------------------------------
# # simple aggregation of the data by cell type
# cts_sample_CellidSmall <- AggregateExpression(object = scobj,
#                                               group.by = c("expertAnno.l1"),
#                                               assays = 'RNA',
#                                               slot = "counts",
#                                               return.seurat = T)
# 
# # add the info to the metadata
# cts_sample_CellidSmall$cell_id <- rownames(cts_sample_CellidSmall@meta.data)
