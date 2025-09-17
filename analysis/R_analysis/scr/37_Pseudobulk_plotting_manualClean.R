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

# DESeq2 processing -------------------------------------------------------
# focus on the immune cells
# 1. Get counts matrix
counts_imm <- cts$RNA %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  dplyr::select("gene",contains("IMMUNE")) %>% 
  column_to_rownames("gene")

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_imm))

colData <- colData %>%
  separate(samples,into = c("tissue","pathology","cell_type"),remove = F,sep = "_")

# get more information from metadata

# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_imm,
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

pdf("../../out/image/heatmap_ht3_shr_IMMMUNE.pdf",width = 20,height = 3) 
draw(ht3_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# select only a subset of genes for the heatmap
# define the genes of interest
GOI <- df_GOI_short %>% 
  filter(Cell_type == "IMMUNE") %>% 
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

pdf("../../out/image/heatmap_ht3_shr_IMMUNE_short.pdf",width = 7,height = 3.5) 
draw(ht3_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# use all samples ---------------------------------------------------------
# focus on the immune cells
# 1. Get counts matrix
counts_all <- cts$RNA %>% 
  data.frame()

# 2. generate sample level metadata
colData_all <- data.frame(samples = colnames(counts_all))

colData_all <- colData_all %>%
  separate(samples,into = c("tissue","pathology","cell_type"),remove = F,sep = "_")

# get more information from metadata

# perform DESeq2 --------
# Create DESeq2 object   
dds_all <- DESeqDataSetFromMatrix(countData = counts_all,
                              colData = colData_all,
                              design = ~ tissue+cell_type+pathology)

# filter
keep_all <- rowSums(counts(dds_all)) >=10
dds_all_filter <- dds_all[keep_all,]

# scale the data
vds_all_filter <- vst(dds_all_filter, blind = F)

# # run DESeq2
# dds_all_filter_tot <- DESeq(dds_all_filter)
# 
# # Check the coefficients for the comparison
# resultsNames(dds_all_filter_tot)
# 
# # Generate results object
# res_all <- results(dds_filter_tot, name = "tissue_WM_vs_CX")
# res_all %>% 
#   data.frame() %>% 
#   rownames_to_column("gene") %>% 
#   arrange(padj)

# plot heatmap ------------------------------------------------------------
# define the genes of interest
gene_id_all <- rownames(assay(vds_all_filter)) %in% GOI

# how many are in the dataset
sum(gene_id_all)

# define the matrices
mat_all <- assay(vds_all_filter)[gene_id_all, ]
mat2_all <- (mat_all - rowMeans(mat_all))/rowSds(mat_all)

#
meta_sample_all <- data.frame(colname = colnames(mat2_all)) %>% 
  left_join(colData_all,by=c("colname"="samples"))

# make the column of the matrix more readable
colnames(mat2_all) <- meta_sample_all$colname

row_ha_all <- rowAnnotation(tissue = meta_sample_all$tissue,
                            pathology = meta_sample_all$pathology,
                            cell_type = meta_sample_all$cell_type,
                            col = list(tissue = c("CX" = "black", "WM" = "gray")))

ht3_all_shr <- Heatmap(t(mat2_all), show_column_names = T,
                       name = "exp",
                       column_title = "senmayo siganture",
                       right_annotation = row_ha_all) 

pdf("../../out/image/heatmap_ht3_shr_ALL.pdf",width = 20,height = 15) 
draw(ht3_all_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# try to plot the comparisons ---------------------------------------------
# plot the comparisons WM CA vs cortex control as FC
# focus on the immune cells
# 1. Get counts matrix
counts_imm_sample <- cts_sample$RNA %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  dplyr::select("gene",contains("IMMUNE")) %>% 
  dplyr::select("gene",contains(c("CX_Ctrl","WM_CA"))) %>% 
  column_to_rownames("gene")

# 2. generate sample level metadata
colData_imm_sample <- data.frame(samples = colnames(counts_imm_sample)) %>%
  separate(samples,into = c("sample","tissue","pathology","cell_type"),remove = F,sep = "_") %>% 
  mutate(pathology = factor(pathology,levels = c("Ctrl","CA")))

# get more information from metadata

# perform DESeq2 --------
# Create DESeq2 object   
dds_imm_sample <- DESeqDataSetFromMatrix(countData = counts_imm_sample,
                              colData = colData_imm_sample,
                              design = ~ pathology)

# filter
keep_imm_sample <- rowSums(counts(dds_imm_sample)) >=10
dds_imm_sample_filter <- dds_imm_sample[keep_imm_sample,]

# scale the data
vds_imm_sample_filter <- vst(dds_imm_sample_filter, blind = F)

# run DESeq2
dds_imm_sample_filter_tot <- DESeq(dds_imm_sample_filter)

# Check the coefficients for the comparison
resultsNames(dds_imm_sample_filter_tot)

# Generate results object
res_imm_sample <- results(dds_imm_sample_filter_tot, name = "pathology_CA_vs_Ctrl")
res_imm_sample %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  arrange(padj)

# plot heatmap ------------------------------------------------------------
# define the genes of interest
gene_id_imm_sample <- rownames(assay(vds_imm_sample_filter)) %in% GOI

# how many are in the dataset
sum(gene_id_imm_sample)

# define the matrices
mat_imm_sample <- assay(vds_imm_sample_filter)[gene_id_imm_sample, ]
mat2_imm_sample <- (mat_imm_sample - rowMeans(mat_imm_sample))/rowSds(mat_imm_sample)

#
meta_imm_sample <- data.frame(colname = colnames(mat2_imm_sample)) %>% 
  left_join(colData_imm_sample,by=c("colname"="samples"))

# make the column of the matrix more readable
colnames(mat2_imm_sample) <- meta_imm_sample$colname

row_ha_imm_sample <- rowAnnotation(tissue = meta_imm_sample$tissue,
                            pathology = meta_imm_sample$pathology,
                            # cell_type = mat2_imm_sample$cell_type,
                            col = list(tissue = c("CX" = "black", "WM" = "gray"),
                                       pathology = c("Ctrl" = "green", "CA" = "red")))

ht3_imm_sample <- Heatmap(t(mat2_imm_sample), show_column_names = T,
                   name = "exp", 
                   column_title = "senmayo siganture",right_annotation = row_ha_imm_sample,
                   # row_names_gp = gpar(fontsize = 3),
                   # cluster_rows = F, 
                   # right_annotation = row_ha, 
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                   
) 

pdf("../../out/image/heatmap_ht3_shr_IMMMUNE_sample.pdf",width = 20,height = 4) 
draw(ht3_imm_sample,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# show only the FC
# define the genes of interest
# define the matrices in this case I am plotting the logFC, therefore there is no need for the normalization
mat_imm_sample2 <- res_imm_sample %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>% 
  filter(gene %in% GOI) %>% 
  dplyr::select(gene,IMMUNE_WMCA_vsCXCtrl=log2FoldChange) %>% 
  column_to_rownames("gene")

ht4_imm_sample <- Heatmap(t(mat_imm_sample2), show_column_names = T,
                          name = "log2FoldChange", 
                          column_title = "senmayo siganture",show_column_dend = F
                          # row_names_gp = gpar(fontsize = 3),
                          # cluster_rows = F, 
                          # right_annotation = row_ha, 
                          # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                          
) 

pdf("../../out/image/heatmap_ht4_shr_IMMMUNE_sample.pdf",width = 20,height = 3) 
draw(ht4_imm_sample,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()

# plot all IMMUNE comparisons ---------------------------------------------
# plot the comparisons WM CA vs cortex control as FC
# focus on the immune cells
# 1. Get counts matrix
counts_imm_sample2 <- cts_sample$RNA %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  dplyr::select("gene",contains("IMMUNE")) %>% 
  column_to_rownames("gene")

# 2. generate sample level metadata
colData_imm_sample2 <- data.frame(samples = colnames(counts_imm_sample2)) %>%
  separate(samples,into = c("sample","tissue","pathology","cell_type"),remove = F,sep = "_") %>% 
  mutate(test = paste0(tissue,"_",pathology)) %>%
  mutate(test = factor(test,levels = c("CX_Ctrl","WM_Core","WM_CI","WM_CA","WM_Ctrl","WM_NAWM","CX_Demye","CX_Mye")))

# get more information from metadata

# perform DESeq2 --------
# Create DESeq2 object   
dds_imm_sample2 <- DESeqDataSetFromMatrix(countData = counts_imm_sample2,
                                         colData = colData_imm_sample2,
                                         design = ~ test)

# filter
keep_imm_sample2 <- rowSums(counts(dds_imm_sample2)) >=10
dds_imm_sample_filter2 <- dds_imm_sample2[keep_imm_sample2,]

# scale the data
vds_imm_sample_filter2 <- vst(dds_imm_sample_filter2, blind = F)

# run DESeq2
dds_imm_sample_filter_tot2 <- DESeq(dds_imm_sample_filter2)

# Check the coefficients for the comparison
resultsNames(dds_imm_sample_filter_tot2)

# Generate results object
tests <- c("test_WM_Core_vs_CX_Ctrl",
           "test_WM_CI_vs_CX_Ctrl",
           "test_WM_CA_vs_CX_Ctrl",
           "test_WM_Ctrl_vs_CX_Ctrl",
           "test_WM_NAWM_vs_CX_Ctrl",
           "test_CX_Demye_vs_CX_Ctrl",
           "test_CX_Mye_vs_CX_Ctrl")

list_res_imm_sample <- map(tests,function(x){
        results(dds_imm_sample_filter_tot2, name = x) %>% 
          data.frame() %>%
          rownames_to_column("gene") %>%
          arrange(padj)
      }) %>% 
  setNames(tests)

#
df_test <- pmap(list(list_res_imm_sample,names(list_res_imm_sample)),function(x,name){
  x %>% 
    dplyr::select("gene","log2FoldChange") %>% 
    dplyr::rename(!!name := "log2FoldChange")
    
}) %>% 
  purrr::reduce(left_join,by="gene")

# show only the FC
# define the genes of interest
# define the matrices in this case I am plotting the logFC, therefore there is no need for the normalization
mat_imm_sample22 <- df_test %>% 
  filter(gene %in% GOI) %>% 
  column_to_rownames("gene")

ht5_imm_sample <- Heatmap(t(mat_imm_sample22), show_column_names = T,
                          name = "log2FoldChange", 
                          column_title = "senmayo siganture"
                          # row_names_gp = gpar(fontsize = 3),
                          # cluster_rows = F, 
                          # right_annotation = row_ha, 
                          # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                          
) 

pdf("../../out/image/heatmap_ht5_shr_IMMMUNE_sample.pdf",width = 20,height = 4) 
draw(ht5_imm_sample,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()
