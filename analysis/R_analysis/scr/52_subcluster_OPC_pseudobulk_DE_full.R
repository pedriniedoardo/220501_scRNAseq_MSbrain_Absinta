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
scobj <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegAllSoupX_expertAnno.rds")

DimPlot(scobj,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(scobj,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(scobj,label = T,raster = T,group.by = "pathology_class")

# wrangling ---------------------------------------------------------------
# explore the full meta
scobj@meta.data
# aggregate the expression per sample per treatment, donor and cell type
cts_sample_all <- AggregateExpression(object = scobj,
                                      group.by = c("orig.ident","pathology_class","origin","expertAnno.l1"),
                                      assays = 'RNA',
                                      slot = "counts",
                                      return.seurat = FALSE)

colnames(cts_sample_all$RNA %>% data.frame())

# OPC processing -----------------------------------------------------------
# focus on the OPC cluster of cells
# 1. Get counts matrix
counts_OPC_Ctrl <- cts_sample_all$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  # pull only the OPC cells
  dplyr::select("gene",contains("OPC")) %>%
  # focus only of the OPC from ctrl samples
  dplyr::select("gene",contains("Ctrl")) %>%
  # remove the unassigned donors
  dplyr::select(-contains("unassigned")) %>%
  column_to_rownames("gene")
#
counts_OPC_Mye <- cts_sample_all$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  # pull only the OPC cells
  dplyr::select("gene",contains("OPC")) %>%
  # focus only of the OPC from ctrl samples
  dplyr::select("gene",contains("CX_Mye")) %>%
  # remove the unassigned donors
  dplyr::select(-contains("unassigned")) %>%
  column_to_rownames("gene")
#
counts_OPC_NAWM <- cts_sample_all$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  # pull only the OPC cells
  dplyr::select("gene",contains("OPC")) %>%
  # focus only of the OPC from ctrl samples
  dplyr::select("gene",contains("NAWM")) %>%
  # remove the unassigned donors
  dplyr::select(-contains("unassigned")) %>%
  column_to_rownames("gene")

# confirm the tables are registrated
identical(rownames(counts_OPC_Ctrl),rownames(counts_OPC_Mye))
identical(rownames(counts_OPC_Mye),rownames(counts_OPC_NAWM))

# merge the three tables
counts_OPC <- bind_cols(counts_OPC_Ctrl,counts_OPC_Mye,counts_OPC_NAWM)

# 2. generate sample level metadata
LUT_sample <- scobj@meta.data %>%
  group_by(orig.ident,sample_id,facility,origin,disease,sex,age,pathology,patient,plaque,PMI,origin_alt) %>%
  summarise()

colData_OPC <- data.frame(samples = colnames(counts_OPC)) %>%
  separate(samples,into = c(c("orig.ident","origin_01","pathology_class","origin_02","expertAnno.l1")),remove = F,sep = "_") %>%
  left_join(LUT_sample,by = c("orig.ident"))

# save matrix and metadata
saveRDS(counts_OPC,file = "../../out/object/counts_OPC_full.rds")
saveRDS(colData_OPC,file = "../../out/object/colData_OPC_full.rds")

# perform DESeq2 ----------------------------------------------------------
# build the model
location_full <- colData_OPC$origin_01
treat_full <- colData_OPC$disease
design_OPC <- model.matrix(~ location_full+treat_full)
colnames(design_OPC)[1] <- c("intercept")
# colnames(design_OPC)[4] <- c("location_fullWM.treat_fullMS")
# save the disign
saveRDS(design_OPC,"../../out/object/design_OPC_ctrl_refCX_full.rds")

# Create DESeq2 object
dds_OPC <- DESeqDataSetFromMatrix(countData = counts_OPC,
                                    colData = colData_OPC,
                                    design = design_OPC)

# filter
# filter at least 10 read per gene per sample
# edgeR::filterByExpr(y1, group = bulk$type)
keep_OPC <- rowSums(counts(dds_OPC)) >=10
# keep_OPC <- rowSums(counts(dds_OPC)) >=10
dds_OPC_filter <- dds_OPC[keep_OPC,]

# plot the raw number of reads per sample
colSums(counts(dds_OPC_filter)) %>%
  data.frame(tot_reads = .) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=sample,y=tot_reads)) + geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))

# scale the data
vds_OPC_filter <- vst(dds_OPC_filter, blind = F)

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/image/heatmap_SampleCluster_pseudobulk_OPC_ctrl_refCX_full.pdf",width = 10,height = 6)
set.seed(1)
hm_OPC <- plot_sample_clustering(vds_OPC_filter,
                                   anno_vars = c("origin_01","disease"),
                                   distance = "euclidean")
draw(hm_OPC,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

# PCA ---------------------------------------------------------------------
plot_vsd_OPC <- plotPCA(vds_OPC_filter,
                          intgroup = c("origin_01","disease","facility")) +
  theme_bw()

plot_vsd_OPC$data %>%
  mutate(sample = str_extract(name,pattern = "s\\d+")) %>%
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=origin_01,label = sample,shape=disease)) +
  geom_point(size =3,alpha=0.6) +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd_OPC$labels[1]) + xlab(plot_vsd_OPC$labels[2])
ggsave("../../out/image/PCA_pseudobulk_OPC_ctrl_refCX_full.pdf",width = 6,height = 4)

plot_vsd_OPC$data %>%
  mutate(sample = str_extract(name,pattern = "s\\d+")) %>%
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=disease,label = sample,shape=origin_01)) +
  geom_point(size =3,alpha=0.6) +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd_OPC$labels[1]) + xlab(plot_vsd_OPC$labels[2])
ggsave("../../out/image/PCA_pseudobulk_OPC_ctrl_refCX_full2.pdf",width = 6,height = 4)

# cell composition --------------------------------------------------------
# martina suggested to check the relative proportion of different cell types from different sample to possibly justify the reason why this sample seems to be missclassified
sample_prop_wide <- scobj@meta.data %>%
  group_by(orig.ident,expertAnno.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  # make it long
  dplyr::select(orig.ident,expertAnno.l1,prop) %>%
  pivot_wider(names_from = orig.ident,values_from = prop,values_fill = 0) %>%
  column_to_rownames("expertAnno.l1") %>%
  # select only the samples in the dataset
  dplyr::select(colData_OPC$orig.ident)

# plot the data as heatmap
meta_sample_prop <- data.frame(colname = colnames(sample_prop_wide)) %>%
  left_join(colData_OPC,by=c("colname"="orig.ident"))

column_meta_sample_prop <- HeatmapAnnotation(location = meta_sample_prop$origin_01,
                                             disease = meta_sample_prop$disease,
                                             col = list(location = c("WM" = "blue",
                                                                     "CX" = "yellow"),
                                                        disease = c("CTRL" = "green",
                                                                    "MS" = "red")))

ht2_shr_OPC2 <- Heatmap(sample_prop_wide, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                          name = "prop_cell_type",
                          column_title = "sample",
                          col = viridis::viridis(option = "turbo",n = 10),
                          
                          # row_names_gp = gpar(fontsize = 3),
                          top_annotation = column_meta_sample_prop,show_row_names = T
                          # cluster_rows = F,
                          # right_annotation = row_ha,
                          # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                          
)
pdf("../../out/image/heatmap_cellType_pseudobulk_OPC_ctrl_refCX_full.pdf",width = 7,height = 5)
draw(ht2_shr_OPC2,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# run DE ------------------------------------------------------------------
# run DESeq2
ddsHTSeq_OPC_filter <- DESeq(dds_OPC_filter)

# Check the coefficients for the comparison
resultsNames(ddsHTSeq_OPC_filter)

# save the filtered object
saveRDS(ddsHTSeq_OPC_filter,"../../out/object/ddsHTSeq_pseudobulk_OPC_ctrl_refCX_full.rds")

# print the contrast
resultsNames(ddsHTSeq_OPC_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast_OPC <- makeContrasts(OPC_Ctrl_WM_vs_CX = location_fullWM,
                                OPC_MS_vs_CTRL = treat_fullMS,
                                # OPC_WM_MS_int = location_fullWM.treat_fullMS,
                                levels = design_OPC)

# pull the results table
res_OPC.WM.Ctrl <- results(ddsHTSeq_OPC_filter, contrast=contrast_OPC[,"OPC_Ctrl_WM_vs_CX"],alpha = 0.05)
res_OPC.MS <- results(ddsHTSeq_OPC_filter, contrast=contrast_OPC[,"OPC_MS_vs_CTRL"],alpha = 0.05)
# res_OPC.WM.MS_int <- results(ddsHTSeq_OPC_filter, contrast=contrast_OPC[,"OPC_WM_MS_int"],alpha = 0.05)
summary(res_OPC.WM.Ctrl)
summary(res_OPC.MS)
# summary(res_OPC.WM.MS_int)

# add the gene symbols
df_res <-
  list(OPC.WM.Ctrl = res_OPC.WM.Ctrl) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCTRL")
# save the table
df_res %>%
  write_tsv("../../out/table/DE_pseudobulk_OPC_ctrl_refCX_cull.tsv")

# add the gene symbols
df_res2 <-
  list(OPC.MS = res_OPC.MS) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCTRL")
# save the table
df_res2 %>%
  write_tsv("../../out/table/DE_pseudobulk_OPC_MS_full.tsv")

# shrink ------------------------------------------------------------------
res_OPC.WM.Ctrl_shr <- lfcShrink(ddsHTSeq_OPC_filter, res = res_OPC.WM.Ctrl, type = "ashr")
res_OPC.MS_shr <- lfcShrink(ddsHTSeq_OPC_filter, res = res_OPC.MS, type = "ashr")

# save the table of DEGs
df_res_shr <-
  list(OPC.WM.Ctrl_shr = res_OPC.WM.Ctrl_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCTRL")
# save the table
df_res_shr %>%
  write_tsv("../../out/table/DE_pseudobulk_OPC_ctrl_refCX_shr_full.tsv")

# save the table of DEGs
df_res_shr2 <-
  list(OPC.MS_shr = res_OPC.MS_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCTRL")
# save the table
df_res_shr2 %>%
  write_tsv("../../out/table/DE_pseudobulk_OPC_MS_shr_full.tsv")


# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
df_res %>%
  data.frame()%>%
  dplyr::filter(baseMean>1) %>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  facet_wrap(~conditionVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../../out/image/histogram_pvalue_pseudobulk_OPC_ctrl_refCX_full.pdf",width = 6,height = 4)

df_res2 %>%
  data.frame()%>%
  dplyr::filter(baseMean>1) %>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  facet_wrap(~conditionVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../../out/image/histogram_pvalue_pseudobulk_OPC_MS_full.pdf",width = 6,height = 4)

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
    data = plot_volcano %>% group_by(conditionVsCTRL) %>% arrange(padj) %>% dplyr::slice(1:100),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/image/vulcano_plot_pseudobulk_OPC_ctrl_refCX_full.pdf",width = 10,height = 10)

# add the info of the genename
plot_volcano2 <- df_res2 %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano2 %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano2[plot_volcano2$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano2[plot_volcano2$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
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
    data = plot_volcano2 %>% group_by(conditionVsCTRL) %>% arrange(padj) %>% dplyr::slice(1:50) %>% filter(abs(log2FoldChange)>1),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/image/vulcano_plot_pseudobulk_OPC_MS_full.pdf",width = 10,height = 10)

# #
# plot_volcano_shr <- df_res_shr %>%
#   # add a clor variable in case significant
#   mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
#   filter(!is.na(col))
# 
# plot_volcano_shr %>%
#   ggplot(aes(x=log2FoldChange,y=-log(padj)))+
#   # geom_point()
#   geom_point(data = plot_volcano_shr[plot_volcano_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
#   geom_point(data = plot_volcano_shr[plot_volcano_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
#   geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
#   geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
#   scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
#   ggrepel::geom_text_repel(
#     data = plot_volcano_shr %>% group_by(conditionVsCX) %>% arrange(padj) %>% dplyr::slice(1:10),
#     aes(label = symbol),segment.alpha=0.4) +
#   facet_wrap(~conditionVsCX)+
#   theme_bw()+
#   theme(strip.background = element_blank())+
#   theme(legend.position = "none")
# ggsave("../../out/image/vulcano_plot_pseudobulk_OPC_ctrl_refCX_shr.pdf",width = 10,height = 10)

# MA plot -----------------------------------------------------------------
# df_res %>%
#   filter(!is.na(padj)) %>%
#   data.frame()%>%
#   mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
#   mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
#   arrange(color) %>%
#   ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
#   scale_x_log10() + scale_color_manual(values = c("gray","red")) +
#   facet_wrap(~conditionVsCTRL)+
#   theme_bw()+
#   theme(strip.background = element_blank()) +
#   theme(legend.position = "none")
# ggsave("../../out/image/MA_plot_pseudobulk_OPC_ctrl_refCX_test.pdf",width = 10,height = 10)

# # shr
# df_res_shr %>%
#   filter(!is.na(padj)) %>%
#   data.frame()%>%
#   mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
#   mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
#   arrange(color) %>%
#   ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
#   scale_x_log10() + scale_color_manual(values = c("gray","red")) +
#   facet_wrap(~conditionVsCTRL)+
#   theme_bw()+
#   theme(strip.background = element_blank()) +
#   theme(legend.position = "none")
# ggsave("../../out/image/MA_plot_pseudobulk_OPC_ctrl_refCX_shr_test.pdf",width = 10,height = 10)

# heatmaps ----------------------------------------------------------------
# pull the scaled values
mat_filter_OPC <- assay(vds_OPC_filter) %>%
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

mat_shr_OPC <- mat_filter_OPC[rownames(vds_OPC_filter) %in% DEG_2, ]
mat2_shr_OPC <- (mat_shr_OPC - rowMeans(mat_shr_OPC))/rowSds(mat_shr_OPC,useNames = TRUE)
#
meta_sample_OPC <- data.frame(colname = colnames(mat2_shr_OPC)) %>%
  left_join(colData_OPC,by=c("colname"="samples"))

# make the column of the matrix more readable
colnames(mat2_shr_OPC) <- meta_sample_OPC$colname

column_ha_shr_OPC <- HeatmapAnnotation(location = meta_sample_OPC$origin_01,
                                         disease = meta_sample_OPC$disease,
                                         col = list(location = c("WM" = "blue",
                                                                 "CX" = "yellow"),
                                                    disease = c("CTRL" = "green",
                                                                "MS" = "red")))

ht2_shr_OPC <- Heatmap(mat2_shr_OPC, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                         name = "exp",
                         column_title = "OPC_shr",
                         # row_names_gp = gpar(fontsize = 3),
                         top_annotation = column_ha_shr_OPC,show_row_names = F,
                         # cluster_rows = F,
                         # right_annotation = row_ha,
                         # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                         
)
pdf("../../out/image/heatmap_DEG_plot_pseudobulk_OPC_ctrl_refCX_shr_full.pdf",width = 7,height = 14)
draw(ht2_shr_OPC,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

#
# the DEGs plot
DEG_3 <- df_res_shr2 %>%
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

mat_shr_OPC2 <- mat_filter_OPC[rownames(vds_OPC_filter) %in% DEG_3, ]
mat2_shr_OPC2 <- (mat_shr_OPC2 - rowMeans(mat_shr_OPC2))/rowSds(mat_shr_OPC2,useNames = TRUE)
#
meta_sample_OPC2 <- data.frame(colname = colnames(mat2_shr_OPC2)) %>%
  left_join(colData_OPC,by=c("colname"="samples"))

# make the column of the matrix more readable
colnames(mat2_shr_OPC2) <- meta_sample_OPC2$colname

column_ha_shr_OPC2 <- HeatmapAnnotation(location = meta_sample_OPC2$origin_01,
                                          disease = meta_sample_OPC2$disease,
                                          col = list(location = c("WM" = "blue",
                                                                  "CX" = "yellow"),
                                                     disease = c("CTRL" = "green",
                                                                 "MS" = "red")))

ht2_shr_OPC2 <- Heatmap(mat2_shr_OPC2, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                          name = "exp",
                          column_title = "OPC_shr",
                          # row_names_gp = gpar(fontsize = 3),
                          top_annotation = column_ha_shr_OPC2,show_row_names = F,
                          # cluster_rows = F,
                          # right_annotation = row_ha,
                          # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                          
)
pdf("../../out/image/heatmap_DEG_plot_pseudobulk_OPC_MS_shr_full.pdf",width = 7,height = 14)
draw(ht2_shr_OPC2,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()


# # upset plot --------------------------------------------------------------
# # read in the table of DEGs
# df_res_shr <- read_tsv("../../out/table/DE_pseudobulk_OPC_ctrl_refCX_shr.tsv")
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
# pdf("../../out/image/upset_DEG_UP_plot_pseudobulk_OPC_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_up), order.by = "freq",nsets = 7)
# dev.off()
# 
# pdf("../../out/image/upset_DEG_DOWN_plot_pseudobulk_OPC_refBASELINE_shr.pdf",width = 14,height = 7)
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
#   write_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_OPC_refBASELINE_shr.tsv")
# 
# df_int_DOWN %>%
#   write_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_OPC_refBASELINE_shr.tsv")
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
pdf("../../out/image/dispersion_plot_pseudobulk_OPC_ctrl_refCX_full.pdf",width = 5,height = 5)
plotDispEsts(ddsHTSeq_OPC_filter)
dev.off()
