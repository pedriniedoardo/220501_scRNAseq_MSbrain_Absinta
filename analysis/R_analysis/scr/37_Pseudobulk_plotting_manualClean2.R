# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(GGally)
library(cowplot)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(DESeq2)

# read in the final object ------------------------------------------------
# scobj <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegration_AllSoupX.rds")
# this one is the dataset using 4000 HVG, which Martina recommended using.
scobj <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegration_AllSoupX_test.rds")

# add the cell_id based on martina's annotation
scobj$cell_id <- scobj@meta.data %>%
  mutate(cell_id = case_when(seurat_clusters %in% c(0,4,14)~"OLIGO",
                             seurat_clusters %in% c(9)~"OPC",
                             seurat_clusters %in% c(3,12)~"ASTRO",
                             seurat_clusters %in% c(5)~"IMMUNE",
                             seurat_clusters %in% c(13)~"LYM",
                             seurat_clusters %in% c(11)~"VAS",
                             seurat_clusters %in% c(1, 2, 10,6)~"EXC NEU",
                             seurat_clusters %in% c(7,8)~"INH NEU")) %>% 
  pull(cell_id)

# filter only the genes in the senmayo signature
list_sig <- readRDS("../../data/signatures/senescence_pathways.rds")
GOI <- list_sig$senmayo %>% 
  pull(Genes)

df_GOI_short <- read_csv("../../data/20231017_Heatmap_SenMayo.csv")

# wrangling ---------------------------------------------------------------
# calculate the average expression per sample per tissue type
cts <- AggregateExpression(object = scobj,
                           group.by = c("pathology_class", "cell_id"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts_sample <- AggregateExpression(object = scobj,
                                  group.by = c("orig.ident","pathology_class", "cell_id"),
                                  assays = 'RNA',
                                  slot = "counts",
                                  return.seurat = FALSE)


# astro processing --------------------------------------------------------
# focus on the immune cells
# 1. Get counts matrix
counts_astro <- cts$RNA %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  dplyr::select("gene",contains("ASTRO")) %>% 
  column_to_rownames("gene")

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_astro))

colData <- colData %>%
  separate(samples,into = c("tissue","pathology","cell_type"),remove = F,sep = "_")

# get more information from metadata

# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_astro,
                              colData = colData,
                              design = ~ tissue)

# filter
keep <- rowSums(counts(dds)) >=10
dds_filter <- dds[keep,]

# scale the data
vds_filter <- vst(dds_filter, blind = F)

# # run DESeq2
# dds_filter_tot <- DESeq(dds_filter)
# 
# # Check the coefficients for the comparison
# resultsNames(dds_filter_tot)
# 
# # Generate results object
# res <- results(dds_filter_tot, name = "tissue_WM_vs_CX")
# res %>% 
#   data.frame() %>% 
#   rownames_to_column("gene") %>% 
#   arrange(padj)

# plot heatmap ------------------------------------------------------------
# define the genes of interest
gene_id <- rownames(assay(vds_filter)) %in% GOI

# how many are in the dataset
sum(gene_id)

# define the matrices
mat <- assay(vds_filter)[gene_id, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat)

#
meta_sample <- data.frame(colname = colnames(mat2)) %>% 
  left_join(colData,by=c("colname"="samples"))

# make the column of the matrix more readable
colnames(mat2) <- meta_sample$colname

sample_ordered_shr <- meta_sample$tissue
column_ha_shr <- HeatmapAnnotation(tissue = sample_ordered_shr,  
                                   col = list(treat = c("CX" = "black", "WM" = "gray"))) 

ht2_shr <- Heatmap(mat2, show_column_names = T,
                   name = "exp", 
                   column_title = "senmayo siganture",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                   
) 

row_ha <- rowAnnotation(tissue = sample_ordered_shr,  
                        col = list(tissue = c("CX" = "black", "WM" = "gray")))

ht3_shr <- Heatmap(t(mat2), show_column_names = T,
                   name = "exp", 
                   column_title = "senmayo siganture",right_annotation = row_ha,
                   # row_names_gp = gpar(fontsize = 3),
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                   
) 

pdf("../../out/image/heatmap_ht3_shr_ASTRO.pdf",width = 20,height = 3) 
draw(ht3_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# select only a subset of genes for the heatmap
# define the genes of interest
GOI <- df_GOI_short %>% 
  filter(Cell_type == "ASTRO") %>% 
  pull(Gene)

gene_id <- rownames(assay(vds_filter)) %in% GOI

# how many are in the dataset
sum(gene_id)

# define the matrices
mat <- assay(vds_filter)[gene_id, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat)

#
meta_sample <- data.frame(colname = colnames(mat2)) %>% 
  left_join(colData,by=c("colname"="samples"))

# make the column of the matrix more readable
colnames(mat2) <- meta_sample$colname

sample_ordered_shr <- meta_sample$tissue
column_ha_shr <- HeatmapAnnotation(tissue = sample_ordered_shr,  
                                   col = list(treat = c("CX" = "black", "WM" = "gray"))) 

ht2_shr <- Heatmap(mat2, show_column_names = T,
                   name = "exp", 
                   column_title = "senmayo siganture",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                   
) 

row_ha <- rowAnnotation(tissue = sample_ordered_shr,  
                        col = list(tissue = c("CX" = "black", "WM" = "gray")))

ht3_shr <- Heatmap(t(mat2), show_column_names = T,
                   name = "exp", 
                   column_title = "senmayo siganture",right_annotation = row_ha,
                   # row_names_gp = gpar(fontsize = 3),
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                   
) 

pdf("../../out/image/heatmap_ht3_shr_ASTRO_short.pdf",width = 7,height = 3.5) 
draw(ht3_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# DESeq2 processing -------------------------------------------------------
# focus on the immune cells
# 1. Get counts matrix
counts_vas <- cts$RNA %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  dplyr::select("gene",contains("VAS")) %>% 
  column_to_rownames("gene")

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_vas))

colData <- colData %>%
  separate(samples,into = c("tissue","pathology","cell_type"),remove = F,sep = "_")

# get more information from metadata

# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_vas,
                              colData = colData,
                              design = ~ tissue)

# filter
keep <- rowSums(counts(dds)) >=10
dds_filter <- dds[keep,]

# scale the data
vds_filter <- vst(dds_filter, blind = F)

# # run DESeq2
# dds_filter_tot <- DESeq(dds_filter)
# 
# # Check the coefficients for the comparison
# resultsNames(dds_filter_tot)
# 
# # Generate results object
# res <- results(dds_filter_tot, name = "tissue_WM_vs_CX")
# res %>% 
#   data.frame() %>% 
#   rownames_to_column("gene") %>% 
#   arrange(padj)

# plot heatmap ------------------------------------------------------------
# define the genes of interest
gene_id <- rownames(assay(vds_filter)) %in% GOI

# how many are in the dataset
sum(gene_id)

# define the matrices
mat <- assay(vds_filter)[gene_id, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat)

#
meta_sample <- data.frame(colname = colnames(mat2)) %>% 
  left_join(colData,by=c("colname"="samples"))

# make the column of the matrix more readable
colnames(mat2) <- meta_sample$colname

sample_ordered_shr <- meta_sample$tissue
column_ha_shr <- HeatmapAnnotation(tissue = sample_ordered_shr,  
                                   col = list(treat = c("CX" = "black", "WM" = "gray"))) 

ht2_shr <- Heatmap(mat2, show_column_names = T,
                   name = "exp", 
                   column_title = "senmayo siganture",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                   
) 

row_ha <- rowAnnotation(tissue = sample_ordered_shr,  
                        col = list(tissue = c("CX" = "black", "WM" = "gray")))

ht3_shr <- Heatmap(t(mat2), show_column_names = T,
                   name = "exp", 
                   column_title = "senmayo siganture",right_annotation = row_ha,
                   # row_names_gp = gpar(fontsize = 3),
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                   
) 

pdf("../../out/image/heatmap_ht3_shr_VAS.pdf",width = 20,height = 3) 
draw(ht3_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# select only a subset of genes for the heatmap
# define the genes of interest
GOI <- df_GOI_short %>% 
  filter(Cell_type == "VAS") %>% 
  pull(Gene)

gene_id <- rownames(assay(vds_filter)) %in% GOI

# how many are in the dataset
sum(gene_id)

# define the matrices
mat <- assay(vds_filter)[gene_id, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat)

#
meta_sample <- data.frame(colname = colnames(mat2)) %>% 
  left_join(colData,by=c("colname"="samples"))

# make the column of the matrix more readable
colnames(mat2) <- meta_sample$colname

sample_ordered_shr <- meta_sample$tissue
column_ha_shr <- HeatmapAnnotation(tissue = sample_ordered_shr,  
                                   col = list(treat = c("CX" = "black", "WM" = "gray"))) 

ht2_shr <- Heatmap(mat2, show_column_names = T,
                   name = "exp", 
                   column_title = "senmayo siganture",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                   
) 

row_ha <- rowAnnotation(tissue = sample_ordered_shr,  
                        col = list(tissue = c("CX" = "black", "WM" = "gray")))

ht3_shr <- Heatmap(t(mat2), show_column_names = T,
                   name = "exp", 
                   column_title = "senmayo siganture",right_annotation = row_ha,
                   # row_names_gp = gpar(fontsize = 3),
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                   
) 

pdf("../../out/image/heatmap_ht3_shr_VAS_short.pdf",width = 10,height = 3.5) 
draw(ht3_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

