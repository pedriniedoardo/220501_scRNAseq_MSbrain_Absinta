# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(GGally)
library(cowplot)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(DESeq2)
library(RNAseqQC)
library(limma)
library(ashr)
library(magick)
library(UpSetR)
# 
# read in the final object ------------------------------------------------
#
scobj <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_MG_harmonySkipIntegration.rds")

DimPlot(scobj,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(scobj,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(scobj,label = T,raster = T,group.by = "pathology_class")

# how many cells are we excluding
scobj@meta.data %>%
  filter(disease == "CTRL") %>%
  group_by(orig.ident,seurat_clusters) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  ggplot(aes(x=orig.ident,y=prop,fill=seurat_clusters)) + geom_col()+theme_cowplot()
  
# subset the cells to remove cluster 4 and 5
scobj_filter <- subset(scobj,subset = seurat_clusters %in% c(0,1,2,3,6))
DimPlot(scobj_filter,label = T,raster = T,group.by = "seurat_clusters")

# wrangling ---------------------------------------------------------------
# explore the full meta
scobj@meta.data

# aggregate the expression per sample per treatment, donor and cell type
cts_sample_all <- AggregateExpression(object = scobj_filter,
                                      group.by = c("orig.ident","pathology_class","origin"),
                                      assays = 'RNA',
                                      slot = "counts",
                                      return.seurat = FALSE)

# MG processing -----------------------------------------------------------
# focus on the MG cluster of cells
# 1. Get counts matrix
counts_MG <- cts_sample_all$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  # pull only the MG cells
  # dplyr::select("gene",contains("MG")) %>%
  # focus only of the MG from ctrl samples
  dplyr::select("gene",contains("Ctrl")) %>%
  # remove the unassigned donors
  dplyr::select(-contains("unassigned")) %>%
  column_to_rownames("gene")

# 2. generate sample level metadata
LUT_sample <- scobj@meta.data %>%
  group_by(orig.ident,sample_id,facility,origin,disease,sex,age,pathology,patient,plaque,PMI,origin_alt) %>%
  summarise()

colData_MG <- data.frame(samples = colnames(counts_MG)) %>%
  separate(samples,into = c(c("orig.ident","origin_01","pathology_class","origin_02")),remove = F,sep = "_") %>%
  left_join(LUT_sample,by = c("orig.ident"))

# save matrix and metadata
saveRDS(counts_MG,file = "../../out/object/counts_MG_clean.rds")
saveRDS(colData_MG,file = "../../out/object/colData_MG_clean.rds")

# perform DESeq2 ----------------------------------------------------------
# build the model
treat_full <- colData_MG$origin_01
design_MG <- model.matrix(~ treat_full)
colnames(design_MG)[1] <- c("intercept")
# save the disign
saveRDS(design_MG,"../../out/object/design_MG_ctrl_refCX_clean.rds")

# Create DESeq2 object
dds_MG <- DESeqDataSetFromMatrix(countData = counts_MG,
                                 colData = colData_MG,
                                 design = design_MG)

# filter
# filter at least 10 read per gene per sample
keep_MG <- rowSums(counts(dds_MG)) >=10
# keep_MG <- rowSums(counts(dds_MG)) >=10
dds_MG_filter <- dds_MG[keep_MG,]

# plot the raw number of reads per sample
colSums(counts(dds_MG_filter)) %>%
  data.frame(tot_reads = .) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=sample,y=tot_reads)) + geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))

# scale the data
vds_MG_filter <- vst(dds_MG_filter, blind = F)

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/image/heatmap_SampleCluster_pseudobulk_MG_ctrl_refCX_clean.pdf",width = 10,height = 6)
set.seed(1)
hm_MG <- plot_sample_clustering(vds_MG_filter,
                                anno_vars = c("origin_01"),
                                distance = "euclidean")
draw(hm_MG,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

# hm_MG <- plot_sample_clustering(vds_MG_filter,
#                                 anno_vars = c("origin_01","facility"),
#                                 distance = "euclidean")
# draw(hm_MG,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 30), "mm"))

# PCA ---------------------------------------------------------------------
plot_vsd_MG <- plotPCA(vds_MG_filter,
                       intgroup = c("origin_01","facility")) +
  theme_bw()

plot_vsd_MG$data %>%
  mutate(sample = str_extract(name,pattern = "s\\d+")) %>%
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=treat_full,label = sample)) +
  geom_point(size =3,alpha=0.6) +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd_MG$labels[1]) + xlab(plot_vsd_MG$labels[2])
ggsave("../../out/image/PCA_pseudobulk_MG_ctrl_refCX_clean.pdf",width = 6,height = 4)

# color also by facility
plot_vsd_MG$data %>%
  mutate(sample = str_extract(name,pattern = "s\\d+")) %>%
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=facility,shape=origin_01,label = sample)) +
  geom_point(size =3,alpha=0.6) +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd_MG$labels[1]) + xlab(plot_vsd_MG$labels[2])
ggsave("../../out/image/PCA_pseudobulk_MG_ctrl_refCX2_clean.pdf",width = 6,height = 4)

# run DE ------------------------------------------------------------------
# run DESeq2
ddsHTSeq_MG_filter <- DESeq(dds_MG_filter)

# Check the coefficients for the comparison
resultsNames(ddsHTSeq_MG_filter)

# save the filtered object
saveRDS(ddsHTSeq_MG_filter,"../../out/object/ddsHTSeq_pseudobulk_MG_ctrl_refCX_clean.rds")

# print the contrast
resultsNames(ddsHTSeq_MG_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast_MG <- makeContrasts(MG_Ctrl_WM_vs_CX = treat_fullWM,
                             levels = design_MG)

# pull the results table
res_MG.WM.Ctrl <- results(ddsHTSeq_MG_filter, contrast=contrast_MG[,"MG_Ctrl_WM_vs_CX"],alpha = 0.05)
summary(res_MG.WM.Ctrl)

# add the gene symbols
df_res <-
  list(MG.WM.Ctrl = res_MG.WM.Ctrl) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCX")
# save the table
df_res %>%
  write_tsv("../../out/table/DE_pseudobulk_MG_ctrl_refCX_clean.tsv")

# shrink ------------------------------------------------------------------
res_MG.WM.Ctrl_shr <- lfcShrink(ddsHTSeq_MG_filter, res = res_MG.WM.Ctrl, type = "ashr")

# save the table of DEGs
df_res_shr <-
  list(MG.WM.Ctrl_shr = res_MG.WM.Ctrl_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCX")
# save the table
df_res_shr %>%
  write_tsv("../../out/table/DE_pseudobulk_MG_ctrl_refCX_shr_clean.tsv")

# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
df_res %>%
  data.frame()%>%
  dplyr::filter(baseMean>1) %>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  facet_wrap(~conditionVsCX)+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../../out/image/histogram_pvalue_pseudobulk_MG_ctrl_refCX_clean.pdf",width = 6,height = 4)

# volcano -----------------------------------------------------------------
# add the info of the genename
plot_volcano <- df_res %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  ggrepel::geom_text_repel(
    data = plot_volcano %>% group_by(conditionVsCX) %>% arrange(padj) %>% dplyr::slice(1:10),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCX)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/image/vulcano_plot_pseudobulk_MG_ctrl_refCX_clean.pdf",width = 10,height = 10)

#
plot_volcano_shr <- df_res_shr %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = plot_volcano_shr %>% group_by(conditionVsCX) %>% arrange(padj) %>% dplyr::slice(1:10),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCX)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/image/vulcano_plot_pseudobulk_MG_ctrl_refCX_shr_clean.pdf",width = 10,height = 10)

# MA plot -----------------------------------------------------------------
df_res %>%
  filter(!is.na(padj)) %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~conditionVsCX)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
ggsave("../../out/image/MA_plot_pseudobulk_MG_ctrl_refCX_clean.pdf",width = 10,height = 10)

# shr
df_res_shr %>%
  filter(!is.na(padj)) %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~conditionVsCX)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
ggsave("../../out/image/MA_plot_pseudobulk_MG_ctrl_refCX_shr_clean.pdf",width = 10,height = 10)

# heatmaps ----------------------------------------------------------------
# pull the scaled values
mat_filter_MG <- assay(vds_MG_filter) %>%
  data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

# the DEGs plot
DEG_2 <- df_res_shr %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1) %>%
  pull(symbol) %>%
  unique()

# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()

mat_shr_MG <- mat_filter_MG[rownames(vds_MG_filter) %in% DEG_2, ]
mat2_shr_MG <- (mat_shr_MG - rowMeans(mat_shr_MG))/rowSds(mat_shr_MG,useNames = TRUE)
#
meta_sample_MG <- data.frame(colname = colnames(mat2_shr_MG)) %>%
  left_join(colData_MG,by=c("colname"="samples"))

# make the column of the matrix more readable
colnames(mat2_shr_MG) <- meta_sample_MG$colname

column_ha_shr_MG <- HeatmapAnnotation(treat = meta_sample_MG$origin_01,
                                      col = list(treat = c("WM" = "blue",
                                                           "CX" = "yellow")))

ht2_shr_MG <- Heatmap(mat2_shr_MG, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                      name = "exp",
                      column_title = "MG_shr",
                      # row_names_gp = gpar(fontsize = 3),
                      top_annotation = column_ha_shr_MG,show_row_names = F,
                      # cluster_rows = F,
                      # right_annotation = row_ha,
                      # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                      
)
pdf("../../out/image/heatmap_DEG_plot_pseudobulk_MG_ctrl_refCX_shr_clean.pdf",width = 7,height = 14)
draw(ht2_shr_MG,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# # upset plot --------------------------------------------------------------
# # read in the table of DEGs
# df_res_shr <- read_tsv("../../out/table/DE_pseudobulk_MG_ctrl_refCX_shr.tsv")
# 
# # build a list of common elements belonging to each set of fegs
# list_DE_up <- df_res_shr %>%
#   split(f = .$conditionVsCX) %>%
#   lapply(function(x){
#     x %>%
#       filter(padj < 0.05,log2FoldChange>1) %>%
#       pull(symbol) %>%
#       unique()
#   })
# 
# glimpse(list_DE_up)
# 
# list_DE_down <- df_res_shr %>%
#   split(f = .$conditionVsCX) %>%
#   lapply(function(x){
#     x %>%
#       filter(padj < 0.05,log2FoldChange<(-1)) %>%
#       pull(symbol) %>%
#       unique()
#   })
# glimpse(list_DE_down)
# 
# # try the upset plot version
# # library(UpSetR)
# pdf("../../out/image/upset_DEG_UP_plot_pseudobulk_MG_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_up), order.by = "freq",nsets = 7)
# dev.off()
# 
# pdf("../../out/image/upset_DEG_DOWN_plot_pseudobulk_MG_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_down), order.by = "freq",nsets = 7)
# dev.off()
# 
# # pull the intersections
# df1_UP <- lapply(list_DE_up,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# df1_DOWN <- lapply(list_DE_down,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# head(df1_UP)
# head(df1_DOWN)
# 
# df2_UP <- data.frame(gene=unique(unlist(list_DE_up)))
# df2_DOWN <- data.frame(gene=unique(unlist(list_DE_down)))
# 
# head(df2_UP)
# head(df2_DOWN)
# 
# df_int_UP <- lapply(df2_UP$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1_UP %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
# 
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# 
# df_int_DOWN <- lapply(df2_DOWN$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1_DOWN %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
# 
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# 
# df_int_UP %>%
#   write_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_MG_refBASELINE_shr.tsv")
# 
# df_int_DOWN %>%
#   write_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_MG_refBASELINE_shr.tsv")
# 
# head(df_int_UP,n=20)
# head(df_int_DOWN,n=20)
# 
# df_int_UP %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))
# 
# df_int_DOWN %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))

# PLOT DISPERSION ---------------------------------------------------------
pdf("../../out/image/dispersion_plot_pseudobulk_MG_ctrl_refCX_clean.pdf",width = 5,height = 5)
plotDispEsts(ddsHTSeq_MG_filter)
dev.off()