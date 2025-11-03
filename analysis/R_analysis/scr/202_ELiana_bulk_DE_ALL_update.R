# AIM ---------------------------------------------------------------------
# perform DE analysis for the dataset

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(limma)
# library(AnnotationDbi) 
# library(AnnotationHub)
library(ComplexHeatmap)
library(ashr)
library(magick)
library(pals)
library(scales)
library(UpSetR)

# read in the data --------------------------------------------------------
ddsHTSeq_filter <- readRDS("../../out/object/202_dds_all_filter_update.rds")
design <- readRDS("../../out/object/201_design_all_update.rds")
vds_filter <- readRDS("../../out/object/vds_all_filter_update.rds")

LUT_sample <- colData(ddsHTSeq_filter) %>% 
  data.frame()

# differential expression analyisis ---------------------------------------
ddsHTSeq_filter <- DESeq(ddsHTSeq_filter)
# if needed is possible to check the distributions of the counts before and after the normalizatoin
# boxplot(log(counts(ddsHTSeq_structure,normalized = T)))
# boxplot(log(counts(ddsHTSeq_structure,normalized = F)))

# save the filtered object
saveRDS(ddsHTSeq_filter,"../../out/object/202_dds_filter_all_UPDATE_DESeq.rds")

# print the contrast
resultsNames(ddsHTSeq_filter)

 # save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast <- makeContrasts(MS_vs_CTRL = diseaseMS,
                          levels = design)

res_MS <- results(ddsHTSeq_filter, contrast=contrast[,"MS_vs_CTRL"],alpha = 0.05)

summary(res_MS)

# add the gene symbols
list_df <- 
  list("MS" = res_MS) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
      # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  })

# save the tables
pmap(list(list_df,names(list_df)),function(x,y){
  name <- paste0("res_",y,".txt")
  # name
  write_tsv(x,file = paste0("../../out/table/",name))
})

# shrink ------------------------------------------------------------------  
res_MS_shr <- lfcShrink(ddsHTSeq_filter, res = res_MS, type = "ashr")

summary(res_MS_shr)

list_df_shr <- 
  list("MS_shr" = res_MS_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
      # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  })

# save the tables
pmap(list(list_df_shr,names(list_df_shr)),function(x,y){
  name <- paste0("res_",y,".txt")
  # name
  write_tsv(x,file = paste0("../../out/table/",name))
})

# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
bind_rows(list_df,.id = "comparison") %>%
  data.frame()%>%
  dplyr::filter(baseMean>1)%>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  theme_bw()+
  facet_wrap(~comparison)+theme(strip.background = element_blank())
ggsave("../../out/image/202_histogram_pvalue.pdf",width = 6,height = 3)

# PLOTTING RESULTS --------------------------------------------------------
# add the info of the genename
test_plot <- bind_rows(list_df,.id = "comparison") %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

test_plot %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = test_plot[test_plot$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = test_plot[test_plot$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = test_plot[test_plot$col==1,] %>% arrange(padj),
    aes(label = symbol),max.overlaps = 5,segment.alpha=0.4,
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  facet_wrap(~comparison)+
  theme_bw() +
  theme(legend.position = "none",strip.background = element_blank())
ggsave("../../out/image/202_vulcano_plot_text.pdf",width = 12,height = 12)
#
# test_plot_fake <- list_df$HYPvsNORM %>%
#   data.frame()%>%
#   # add a clor variable in case significant
#   mutate(col=ifelse(((pvalue<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0))
# 
# test_plot_fake %>%
#   ggplot(aes(x=log2FoldChange,y=-log(pvalue)))+
#   # geom_point()
#   geom_point(data = test_plot_fake[test_plot_fake$col==0,],aes(x=log2FoldChange,y=-log(pvalue),col=factor(col)),alpha=0.05)+
#   geom_point(data = test_plot_fake[test_plot_fake$col==1,],aes(x=log2FoldChange,y=-log(pvalue),col=factor(col)),alpha=0.5)+
#   geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
#   geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
#   scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
#   ggrepel::geom_text_repel(
#     data = test_plot_fake[test_plot_fake$col==1,],
#     aes(label = symbol),
#     size = 5,
#     box.padding = unit(0.35, "lines"),
#     point.padding = unit(0.3, "lines")) +
#   theme_bw() +
#   theme(legend.position = "none")
# ggsave("../../out/image/vulcano_plot_text_DIAvsCTRL_fake.pdf",width = 12,height = 12)

#
test_plot_shr <- bind_rows(list_df_shr,.id = "comparison") %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

test_plot_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = test_plot_shr[test_plot_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = test_plot_shr[test_plot_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = test_plot_shr[test_plot$col==1,] %>% arrange(padj),
    aes(label = symbol),max.overlaps = 5,segment.alpha=0.4,
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  theme_bw() +
  facet_wrap(~comparison)+
  theme(legend.position = "none",strip.background = element_blank())
ggsave("../../out/image/202_vulcano_plot_text_shr.pdf",width = 12,height = 12)

# plotMA(res, ylim = c(-5, 5))
bind_rows(list_df,.id = "comparison") %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>% 
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) + 
  scale_x_log10() + scale_color_manual(values = c("gray","red")) + theme_bw() + 
  theme(legend.position = "none",strip.background = element_blank())+facet_wrap(~comparison)
ggsave("../../out/image/202_MA_plot.pdf",width = 10,height = 3)

# using the shrinked values
bind_rows(list_df_shr,.id = "comparison") %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>% 
  ggplot(aes(baseMean,y = log2FoldChange)) + 
  geom_point(aes(col=color),alpha=0.2) + 
  ggrepel::geom_text_repel(
    data = test_plot_shr[test_plot_shr$col==1,],
    aes(label = symbol),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  scale_x_log10() + 
  scale_color_manual(values = c("gray","red")) + 
  theme_bw() +
  theme(legend.position = "none",strip.background = element_blank())+facet_wrap(~comparison)
ggsave("../../out/image/202_MA_plot_shr.pdf",width = 24,height = 12)

# # the DEGs plot stringent
# DEG_1 <- list_df$ %>%
#   data.frame()%>%
#   # add a clor variable in case significant
#   mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
#   filter(col==1) %>%
#   pull(symbol)

mat_filter <- assay(vds_filter) %>%
  as.data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

# mat <- mat_filter[rownames(vds_filter) %in% DEG_1, ]
# mat2 <- (mat - rowMeans(mat))/rowSds(mat)
# #
# 
# sample_ordered <- str_extract(colnames(mat2),pattern = "BMP9|mock")
# column_ha <- HeatmapAnnotation(treat = sample_ordered,  
#                                col = list(treat = c("mock" = "green", "BMP9" = "gray"))) 
# 
# ht2 <- Heatmap(mat2, 
#                name = "exp", 
#                column_title = "BMP9",
#                row_names_gp = gpar(fontsize = 3),
#                top_annotation = column_ha, 
#                # cluster_rows = F, 
#                # right_annotation = row_ha, 
#                # row_split = rep(c(1,2,3,4),c(2,3,4,7))
# ) 
# pdf("out/image/heatmap_BMP9_DEG.pdf",width = 4,height = 15) 
# draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
# dev.off()

# plot the DEGs for Veh20Dynes --------------------------------------------
DEG_MS <- list_df_shr$MS_shr %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1) %>%
  pull(symbol)

mat_shr <- mat_filter[rownames(vds_filter) %in% DEG_MS, ]
# pull only the samples of interested based on the comparison
id_sample <- str_detect(vds_filter@colData$Group,pattern = "CTRL|MS")
mat1_shr <- mat_shr[,id_sample]
mat2_shr <- (mat1_shr - rowMeans(mat1_shr))/rowSds(mat1_shr,useNames = TRUE)
#
meta_sample <- data.frame(colname = colnames(mat2_shr)) %>% 
  left_join(LUT_sample,by=c("colname"="HUGE.ID"))

# make the column of the matrix more readable
# colnames(mat2_shr) <- meta_sample$Sample_ID2

# define a LUT for the clones
# color_id <- alphabet(length(unique(meta_sample$clone)))
# check the colors
# show_col(color_id)

# build the named vector
# names(color_id) <- unique(meta_sample$clone)

column_ha_shr <- HeatmapAnnotation(# clone = meta_sample$clone,
                                   condition = meta_sample$Group,
                                   Treatment = meta_sample$Treatment,
                                   Gender = meta_sample$sex,
                                   col = list(# clone = color_id,
                                              condition = c("MS" = "black", "CTRL" = "gray"),
                                              Treatment = c("ctrl" = "green", "Naive" = "red","TRT" = "orange"),
                                              Gender = c("M" = "blue", "F" = "pink")))  

ht2_shr <- Heatmap(mat2_shr, show_column_names = T,
                   name = "exp", 
                   column_title = "DEG_MS_shr",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr,show_row_names = T,
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                   
) 
pdf("../../out/image/202_heatmap_DEG_shr_MS.pdf",width = 7,height = 8) 
draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()


# PLOT DISPERSION ---------------------------------------------------------
pdf("../../out/image/202_ddsHTSeq_filter_dispersion.pdf",width = 5,height = 5) 
plotDispEsts(ddsHTSeq_filter)
dev.off()
