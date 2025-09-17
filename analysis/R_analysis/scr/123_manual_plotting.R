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
scobj <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")

# filter only the genes in the senmayo signature
list_sig <- readRDS("../../data/signatures/senescence_pathways.rds")
GOI <- list_sig$senmayo %>%
  pull(Genes)

# df_GOI_short <- read_csv("../../data/20231017_Heatmap_SenMayo.csv")


# plotting ----------------------------------------------------------------
# plot the UMAP with the hihglight only for the IMM cells
# DimPlot(scobj,group.by = "expertAnno.l1",raster = T)
df_all <- scobj@meta.data %>%
  data.frame() %>%
  rownames_to_column() %>%
  filter(! expertAnno.l1 %in% c("IMM","LYM")) %>%
  left_join(scobj@reductions$umap@cell.embeddings %>%
              data.frame() %>%
              rownames_to_column(),by = "rowname")

df_IMM <- scobj@meta.data %>%
  data.frame() %>%
  rownames_to_column() %>%
  filter(expertAnno.l1 %in% c("IMM","LYM")) %>%
  left_join(scobj@reductions$umap@cell.embeddings %>%
              data.frame() %>%
              rownames_to_column(),by = "rowname")

df_all %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2)) +
  geom_point(size = 0.1,col="gray") +
  geom_point(data = df_IMM,size = 0.1,col="#60b1bd") +
  theme_cowplot()
ggsave("../../out/image/123_UMAP_highlight_IMM.pdf",width = 4,height = 4)

# highlight the vasc
df_all2 <- scobj@meta.data %>%
  data.frame() %>%
  rownames_to_column() %>%
  filter(! expertAnno.l1 %in% c("VAS")) %>%
  left_join(scobj@reductions$umap@cell.embeddings %>%
              data.frame() %>%
              rownames_to_column(),by = "rowname")

df_VAS <- scobj@meta.data %>%
  data.frame() %>%
  rownames_to_column() %>%
  filter(expertAnno.l1 %in% c("VAS")) %>%
  left_join(scobj@reductions$umap@cell.embeddings %>%
              data.frame() %>%
              rownames_to_column(),by = "rowname")

df_all2 %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2)) +
  geom_point(size = 0.1,col="gray") +
  geom_point(data = df_VAS,size = 0.1,col="#60b1bd") +
  theme_cowplot()
ggsave("../../out/image/123_UMAP_highlight_VAS.pdf",width = 4,height = 4)


# highlight the astro
df_all3 <- scobj@meta.data %>%
  data.frame() %>%
  rownames_to_column() %>%
  filter(! expertAnno.l1 %in% c("AST")) %>%
  left_join(scobj@reductions$umap@cell.embeddings %>%
              data.frame() %>%
              rownames_to_column(),by = "rowname")

df_AST <- scobj@meta.data %>%
  data.frame() %>%
  rownames_to_column() %>%
  filter(expertAnno.l1 %in% c("AST")) %>%
  left_join(scobj@reductions$umap@cell.embeddings %>%
              data.frame() %>%
              rownames_to_column(),by = "rowname")

df_all3 %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2)) +
  geom_point(size = 0.1,col="gray") +
  geom_point(data = df_AST,size = 0.1,col="#60b1bd") +
  theme_cowplot()
ggsave("../../out/image/123_UMAP_highlight_AST.pdf",width = 4,height = 4)

# wrangling ---------------------------------------------------------------
# calculate the average expression per sample per tissue type
cts <- AggregateExpression(object = scobj,
                           group.by = c("pathology_class","expertAnno.l1"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

# DESeq2 processing -------------------------------------------------------
# focus on the immune cells
# 1. Get counts matrix
counts_imm <- cts$RNA %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  dplyr::select("gene",contains("IMM")) %>% 
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
keep <- edgeR::filterByExpr(counts(dds), group = colData$pathology)
# keep_ASTRO <- rowSums(counts(dds_ASTRO)) >=10
dds_filter <- dds[keep,]
# keep <- rowSums(counts(dds)) >=10
# dds_filter <- dds[keep,]

# scale the data
vds_filter <- vst(dds_filter, blind = F)

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

pdf("../../out/image/revision/123_heatmap_ht3_shr_IMM.pdf",width = 20,height = 3) 
draw(ht3_shr,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()
