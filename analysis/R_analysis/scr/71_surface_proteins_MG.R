# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
# read in the table for the deg from cluster wise comparison
df <- read_tsv("../../out/table/DE_pseudobulk_MG_ctrl_refCX_shr.tsv")

# reference: http://wlab.ethz.ch/cspa/#abstract
ref <- readxl::read_excel("../../data/S2_File.xlsx",sheet = 1)

# how many are also in Uniprot
table(ref$`CSPA category`,ref$`UniProt Cell surface`,useNA="ifany")

# how many are CDs
table(ref$`CSPA category`,ref$CD!="no")


# -------------------------------------------------------------------------
# how many are the proteins that are differentially expressed
df_tot <- left_join(df,ref,by=c("symbol"="ENTREZ gene symbol")) %>% 
  filter(!is.na(organism))

df_tot %>% 
  write_tsv("../../out/table/surfaceProt_pseudobulk_MG_ctrl_refCX_shr.tsv")

# plot the expression data to focus only on the surface protein genes

# heatmaps ----------------------------------------------------------------
# save the filtered object
ddsHTSeq_MG_filter <- readRDS("../../out/object/ddsHTSeq_pseudobulk_MG_ctrl_refCX.rds")
colData_MG <- readRDS(file = "../../out/object/colData_MG.rds")

# scale the data
vds_MG_filter <- vst(ddsHTSeq_MG_filter, blind = F)

# plot_vsd_MG <- plotPCA(vds_MG_filter,
#                        intgroup = c("origin_01","facility")) +
#   theme_bw()
# 
# plot_vsd_MG$data %>%
#   mutate(sample = str_extract(name,pattern = "s\\d+")) %>%
#   # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
#   ggplot(aes(x=PC1,y=PC2,col=origin_01,label = sample)) +
#   geom_point(size =3,alpha=0.6) +
#   ggrepel::geom_text_repel(show.legend = F)+
#   scale_x_continuous(expand = expansion(mult = 0.1))+
#   theme_bw() + ylab(plot_vsd_MG$labels[1]) + xlab(plot_vsd_MG$labels[2])

# pull the scaled values
mat_filter_MG <- assay(vds_MG_filter) %>%
  data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

# the DEGs plot
DEG_2 <- df_tot %>%
  data.frame() %>%
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

# column annotation
column_ha_shr_MG <- HeatmapAnnotation(location = meta_sample_MG$origin_01,
                                      gender = meta_sample_MG$sex,
                                      col = list(location = c("WM" = "blue",
                                                             "CX" = "yellow"),
                                                 gender = c("M" = "cyan",
                                                            "F" = "pink")))

# row annotation
LUT_genes <- df_tot %>%
  data.frame() %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1)
  
# pull the order from the matrix
LUT_gene_order <- data.frame(symbol = rownames(mat2_shr_MG)) %>%
  left_join(LUT_genes,by="symbol")

row_ha <- rowAnnotation(confidence = LUT_gene_order$CSPA.category,
                        col = list(confidence = c("1 - high confidence" = "green", "2 - putative" = "yellow","3 - unspecific" = "red")))

ht2_shr_MG <- Heatmap(mat2_shr_MG, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                      name = "exp",
                      column_title = "MG_shr",
                      right_annotation = row_ha, 
                      # row_names_gp = gpar(fontsize = 3),
                      top_annotation = column_ha_shr_MG,show_row_names = T,
                      # cluster_rows = F,
                      # right_annotation = row_ha,
                      # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                      
)
pdf("../../out/image/surfaceProt_heatmap_DEG_plot_pseudobulk_MG_ctrl_refCX.pdf",width = 7,height = 30)
draw(ht2_shr_MG,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# volcano -----------------------------------------------------------------
# add the info of the genename
plot_volcano <- df_tot %>%
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
    # data = plot_volcano %>% group_by(conditionVsCX) %>% arrange(padj) %>% dplyr::slice(1:20) %>% filter(abs(log2FoldChange)>1),
    data = plot_volcano %>% group_by(conditionVsCX) %>% arrange(padj) %>% filter(abs(log2FoldChange)>1) %>% filter(padj<0.05),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCX)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/image/surfaceProt_vulcano_plot_pseudobulk_MG_ctrl_refCX.pdf",width = 10,height = 10)
