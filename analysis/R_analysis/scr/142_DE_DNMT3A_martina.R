# AIM ---------------------------------------------------------------------
# test DE of cluster 8 vs cluster 0 at the sc level and pseudobulk.
# check if the expression of DNMT3A is significanlty disregulated

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(tidyverse)
library(harmony)
library(ggExtra)
library(ComplexUpset)
library(cowplot)
library(UpSetR)
library(DESeq2)

# read in the dataset -----------------------------------------------------
sobj <- readRDS("../../data/all20_immune.rds")
DimPlot(sobj,label = T,group.by = "RNA_snn_res.0.5")

# check the object version
class(sobj@assays$RNA)

# test DGE at single cell -------------------------------------------------
sobj@meta.data %>% pull(RNA_snn_res.0.5) %>% table()
# set the ident
Idents(sobj) <- "RNA_snn_res.0.5"

# run the DGE over the same cell tow cell type. The units are the single cells.
# log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
mono.de <- FindMarkers(sobj, ident.1 = "8", ident.2 = "0", verbose = FALSE,logfc.threshold = 0)

head(mono.de)
# The p-values obtained from this analysis should be interpreted with caution, because these tests treat each cell as an independent replicate and ignore inherent correlations between cells originating from the same sample. Such analyses have been shown to find a large number of false positive associations, as has been demonstrated by Squair et al., 2021, Zimmerman et al., 2021, Junttila et al., 2022, and others. Below, we show how pseudobulking can be used to account for such within-sample correlation.

# check the position of DNMT3A in the comparison
mono.de %>%
  data.frame() %>%
  rownames_to_column("symbol") %>%
  # filter(str_detect(gene,pattern = "DNMT"))
  filter(symbol %in% c("DNMT3A"))

# add the info of the genename
plot_volcano <- mono.de %>%
  data.frame() %>%
  rownames_to_column("symbol") %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((p_val_adj<0.05)&abs(avg_log2FC)>0.5&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano %>%
  ggplot(aes(x=avg_log2FC,y=-log(p_val_adj)))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=avg_log2FC,y=-log(p_val_adj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=avg_log2FC,y=-log(p_val_adj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-0.5,0.5),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # ggrepel::geom_text_repel(
  #   data = plot_volcano %>% 
  #     group_by(condVsCTRL) %>% 
  #     arrange(padj) %>% 
  #     dplyr::slice(1:30) %>% 
  #     dplyr::filter(abs(log2FoldChange)>1) %>% 
  #     dplyr::filter(padj<0.05),
  #   aes(label = symbol),segment.alpha=0.4) +
  # facet_wrap(~condVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
# ggsave("../out/plot/02_vulcano_plot_pseudobulk_filterExp.pdf",width = 10,height = 10)

# plot violin for the expression
VlnPlot(sobj,features = "DNMT3A")

# save the table
saveRDS(mono.de,"../../out/object/142_mono.de.rds")

# test DGE at pseudobulk --------------------------------------------------
# pseudobulk the counts based on donor-condition-celltype
pseudo_sobj <- AggregateExpression(sobj, assays = "RNA", return.seurat = T, group.by = c("sample", "RNA_snn_res.0.5"),slot = "counts")

# confirm the count slot is correctly populated, notice that the slot contains integers
pseudo_sobj@assays$RNA@counts
pseudo_sobj@assays$RNA@data

# add the covariate for the stimulation per cell type
# pseudo_sobj$sample.disease <- paste(pseudo_sobj$dataset, pseudo_sobj$disease, sep = "_")
pseudo_sobj$sample.cell <- 
  pseudo_sobj@meta.data %>%
  rownames_to_column("sample.cell") %>%
  pull("sample.cell")

pseudo_sobj$cell <- 
  pseudo_sobj@meta.data %>%
  separate(sample.cell,into = c("sample","cell"),sep = "_") %>%
  pull(cell)

# explore the dimensionality of the new dataset
pseudo_sobj@meta.data %>%
  group_by(cell) %>%
  summarise(n = n())

# run the DGE over the same cell type for MS vs CTRL This time the unist are the pseudobulks per donor/celltype/stimulus
Idents(pseudo_sobj) <- "cell"
bulk.mono.de <- FindMarkers(object = pseudo_sobj, 
                            ident.1 = "8", 
                            ident.2 = "0",
                            test.use = "DESeq2")
head(bulk.mono.de)

# check the position of DNMT3A in the comparison
bulk.mono.de %>%
  data.frame() %>%
  rownames_to_column("symbol") %>%
  # filter(str_detect(gene,pattern = "DNMT"))
  filter(symbol %in% c("DNMT3A"))

# add the info of the genename
plot_volcano_bulk <- bulk.mono.de %>%
  data.frame() %>%
  rownames_to_column("symbol") %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((p_val_adj<0.05)&abs(avg_log2FC)>0.5&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano_bulk %>%
  ggplot(aes(x=avg_log2FC,y=-log(p_val_adj)))+
  # geom_point()
  geom_point(data = plot_volcano_bulk[plot_volcano_bulk$col==0,],aes(x=avg_log2FC,y=-log(p_val_adj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano_bulk[plot_volcano_bulk$col==1,],aes(x=avg_log2FC,y=-log(p_val_adj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-0.5,0.5),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # ggrepel::geom_text_repel(
  #   data = plot_volcano %>% 
  #     group_by(condVsCTRL) %>% 
  #     arrange(padj) %>% 
  #     dplyr::slice(1:30) %>% 
  #     dplyr::filter(abs(log2FoldChange)>1) %>% 
  #     dplyr::filter(padj<0.05),
  #   aes(label = symbol),segment.alpha=0.4) +
  # facet_wrap(~condVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
# ggsave("../out/plot/02_vulcano_plot_pseudobulk_filterExp.pdf",width = 10,height = 10)

# save the table
saveRDS(bulk.mono.de,"../../out/object/142_bulk.mono.de.rds")

# test DGE pseudobulk deseq2 ----------------------------------------------


















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

# wrangling ---------------------------------------------------------------
# explore the full meta
sobj@meta.data %>%
  group_by(sample,disease) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = disease,values_from = n)

# aggregate the expression per sample per donor.id and stimulation
cts_sample_all <- AggregateExpression(object = sobj,
                                      group.by = c("sample","disease","RNA_snn_res.0.5"),
                                      assays = 'RNA',
                                      slot = "counts",
                                      return.seurat = FALSE)

# processing --------------------------------------------------------------
# extract only the condition of interest homeo-MG|MIMS-iron
# 1. Get counts matrix
counts <- cts_sample_all$RNA %>%
  as.data.frame() %>%
  select(contains("_8")|contains("_0"))

# 2. generate sample level metadata
LUT_sample <- sobj@meta.data %>%
  group_by(sample,disease,RNA_snn_res.0.5) %>%
  summarise() %>%
  # harmonyze the naming of the donor_id
  mutate(sample_cell = paste(sample,disease,RNA_snn_res.0.5,sep = "_")) %>%
  filter(RNA_snn_res.0.5 %in% c("0","8")) %>%
  ungroup()

# match the LUT with the expression matrix
colData <- data.frame(sample.id = colnames(counts)) %>%
  # separate(samples,into = c(c("orig.ident","origin_01","pathology_class","origin_02")),remove = F,sep = "_") %>%
  left_join(LUT_sample,by = c("sample.id" = "sample_cell"))

# perform DESeq2 ----------------------------------------------------------
# build the model
treat <- colData$RNA_snn_res.0.5 %>% str_replace(pattern = "-",replacement = ".")
block_donor <- colData$sample

design <- model.matrix(~ block_donor + treat)
# design <- model.matrix(~ treat)
colnames(design)[1] <- c("intercept")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = design)

# filter low aboundant features
feat_keep <- edgeR::filterByExpr(counts(dds), group = colData$RNA_snn_res.0.5)
dds_filter <- dds[feat_keep,]

# plot the raw number of reads per sample
colSums(counts(dds_filter)) %>%
  data.frame(tot_reads = .) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=sample,y=tot_reads)) + geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))

# scale the data
vds_filter <- vst(dds_filter, blind = F)

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/image/142_heatmap_SampleCluster_filterExp.pdf",width = 10,height = 6)
set.seed(1)
hm <- plot_sample_clustering(vds_filter,
                             anno_vars = c("disease","RNA_snn_res.0.5"),
                             distance = "euclidean")
draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

# PCA ---------------------------------------------------------------------
plot_vsd <- plotPCA(vds_filter,
                    intgroup = c("RNA_snn_res.0.5","disease")) +
  theme_bw()

plot_vsd$data %>%
  ggplot(aes(x=PC1,y=PC2,col=RNA_snn_res.0.5)) +
  geom_point(size =3,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave("../../out/image/142_PCA_pseudobulk_filterExp.pdf",width = 6,height = 4)

# run DE ------------------------------------------------------------------
# run DESeq2
ddsHTSeq_filter <- DESeq(dds_filter)

# Check the coefficients for the comparison
resultsNames(ddsHTSeq_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast <- makeContrasts(clu8_vs_clu0 = treat8,
                          levels = design)

# pull the results table
res_MS <- results(ddsHTSeq_filter, contrast=contrast[,"clu8_vs_clu0"],alpha = 0.05)
summary(res_MS)

# add the gene symbols
df_res <-
  list(clu8 = res_MS) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "condVsclu0")

# save the table
df_res %>%
  write_tsv("../../out/table/142_DE_pseudobulk_filterExp.tsv")

# check the position of DNMT3A in the comparison
df_res %>%
  # filter(str_detect(gene,pattern = "DNMT"))
  filter(symbol %in% c("DNMT3A"))

# shrink ------------------------------------------------------------------
res_MS_shr <- lfcShrink(ddsHTSeq_filter, res = res_MS, type = "ashr")

# save the table of DEGs
df_res_shr <-
  list(clu8_shr = res_MS_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "condVsclu0")

# save the table
df_res_shr %>%
  write_tsv("../../out/table/142_DE_pseudobulk_filterExp_shr.tsv")

# check the position of DNMT3A in the comparison
df_res_shr %>%
  # filter(str_detect(gene,pattern = "DNMT"))
  filter(symbol %in% c("DNMT3A"))

plot_volcano <- df_res %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>0.5&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-0.5,0.5),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # ggrepel::geom_text_repel(
  #   data = plot_volcano %>% 
  #     group_by(condVsCTRL) %>% 
  #     arrange(padj) %>% 
  #     dplyr::slice(1:30) %>% 
  #     dplyr::filter(abs(log2FoldChange)>1) %>% 
  #     dplyr::filter(padj<0.05),
  #   aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~condVsclu0)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
# ggsave("../../out/image/142_vulcano_plot_pseudobulk_filterExp.pdf",width = 10,height = 10)

#
plot_volcano_shr <- df_res_shr %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>0.5&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-0.5,0.5),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano_shr %>% group_by(condVsCTRL) %>% arrange(padj) %>% dplyr::slice(1:10)%>% filter(abs(log2FoldChange)>1),
  #   aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~condVsclu0)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
# ggsave("../../out/image/02_vulcano_plot_pseudobulk_filterExp_shr.pdf",width = 10,height = 10)

# plot the expression value
data <- ddsHTSeq_filter

lut <- colData(data) %>%
  data.frame() %>%
  mutate(sample = paste0("X",sample.id))

# GOI <- c(sig_B$symbol,sig_D$symbol)
GOI <- c("DNMT3A")
# GOI <- subset_genes

MR <- counts(data,normalized=T) %>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  pivot_longer(names_to = "sample",values_to = "exp",-symbol)%>%
  group_by(sample)%>%
  summarise(MR = sum(exp)/10^6)

# plot the data following the methods implemented in the plotCounts funciton from DESeq2
# Normalized counts plus a pseudocount of 0.5 are shown by default.
counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = "sample") %>%
  mutate(count_norm_adj = count + 0.5)%>%
  mutate(treat = fct_relevel(RNA_snn_res.0.5,"0")) %>%
  # mutate(symbol = factor(symbol,levels = c("PROCR","GDF9", "GDF11","TGFB1", "TGFB2", "TGFB3","INHBA", "INHBB","MSTN"))) %>%
  ggplot(aes(x=treat,y = count_norm_adj))+
  geom_boxplot(outlier.shape = NA,linewidth = 1)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6,size = 3)+
  facet_wrap(~symbol,scales = "free",ncol=3)+
  scale_y_continuous(trans = "log1p") +
  # scale_y_log10()+
  theme_bw(base_size = 22,base_rect_size = 2)+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = "sample") %>%
  mutate(count_norm_adj = count + 0.5)%>%
  mutate(treat = fct_relevel(RNA_snn_res.0.5,"0")) %>%
  # mutate(symbol = factor(symbol,levels = c("PROCR","GDF9", "GDF11","TGFB1", "TGFB2", "TGFB3","INHBA", "INHBB","MSTN"))) %>%
  ggplot(aes(x=treat,y = count_norm_adj))+
  geom_boxplot(outlier.shape = NA,linewidth = 1)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6,size = 3)+
  facet_grid(disease~symbol,scales = "free")+
  scale_y_continuous(trans = "log1p") +
  # scale_y_log10()+
  theme_bw(base_size = 22,base_rect_size = 2)+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

