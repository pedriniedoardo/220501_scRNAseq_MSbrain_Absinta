# AIM ---------------------------------------------------------------------
# Task for the rebuttal of the paper from Aletta
# check how the curated FTH1/GPNMB gene panel (lysosomal, iron and lipid handling modules) correlates with FTH1 and GPNMB across MG cells, pooling all samples regardless of source (disease, pathology, patient).
# Expression is averaged per sample, correlated with Spearman, tested for significance (BH-adjusted), and summarized as a heatmap (annotated by gene module) and scatter plots highlighting the significantly correlated genes.

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)

# read in the dataset -----------------------------------------------------
# read in the full dataset from WM MG only (Martina's Nature dataset)
data.combined <- readRDS("../../data/all20_immune.rds")
p_umap_all <- DimPlot(data.combined,label = T,raster = F,group.by = "seurat_clusters")
ggsave(plot = p_umap_all,filename = "../../out/image/custom/05_UMAP_all_clusters.pdf",width = 6,height = 5)

# load the table compiled by Aletta
df_GOI <- read_csv("../../data/20260803_Geneset_Aletta_FTH1_GPNMB.csv") %>%
  # LIPA has two annotation, keep only one
  mutate(rep = duplicated(Gene)) %>%
  filter(rep == FALSE) %>%
  dplyr::select(-rep)

# define the GOI mentioned by Aletta
GOI <- c("FTH1","GPNMB")

# wrangling ---------------------------------------------------------------
# Aletta is interested in only the MG subset, therefore remove the clusters that are not MG
data.combined$is_MG <- !(data.combined$seurat_clusters %in% c(10,9,7,4))
data.combined_filter <- subset(data.combined, subset = is_MG == T)

p_umap_MG <- DimPlot(data.combined_filter,label = T,raster = F,group.by = "seurat_clusters")
ggsave(plot = p_umap_MG,filename = "../../out/image/custom/05_UMAP_MG_filtered.pdf",width = 6,height = 5)

# aggregate the expression for all the genes per sample

# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined_filter$group <- paste0(data.combined_filter$sample,"|",
                                     data.combined_filter$disease,"|",
                                     data.combined_filter$pathology,"|",
                                     data.combined_filter$patient,"|")

# define the ident for the grouping
Idents(data.combined_filter) <- "group"
DefaultAssay(data.combined_filter) <- "RNA"
average_RNA <- AverageExpression(data.combined_filter,group.by = c("group")) %>%
  .$RNA

# save also how many cells are contribution to the average
df_group_n <- data.combined_filter@meta.data %>%
  group_by(group) %>%
  summarise(n = n())
write_tsv(df_group_n,"../../out/table/custom/05_Aletta_correlation_group_cellCounts.tsv")

# run the correlation
average_GOI_RNA <- average_RNA[rownames(average_RNA) %in% c(df_GOI$Gene,GOI),] %>%
  as.data.frame()

# remove potentailly problematic genes
average_GOI_RNA_filter <- average_GOI_RNA %>%
  filter(rowMeans(average_GOI_RNA)>0)

# make a gene vs gene correlation matrix regardless of the sample
corr_mat <- cor(t(average_GOI_RNA_filter),method = "spearman")

# annotate the genes by their module (the two GOI genes have no module of their own)
df_module_ann <- tibble(Gene = rownames(corr_mat)) %>%
  left_join(df_GOI,by = "Gene") %>%
  mutate(Module = ifelse(Gene %in% GOI,"GOI",Module))

id_module <- unique(str_subset(df_module_ann$Module,pattern = "GOI",negate = T))
module_palette <- setNames(pals::trubetskoy(n = length(id_module)),nm = id_module)
# add the GOI category
module_palette["GOI"] <- "black"

# show_col(module_palette)

row_ha <- rowAnnotation(Module = df_module_ann$Module,
                        col = list(Module = module_palette))

# highlight the two GOI labels in red
label_col <- ifelse(rownames(corr_mat) %in% GOI,"red","black")

# full correlation matrix
ht01 <- Heatmap(corr_mat,
                left_annotation = row_ha,
                name = "Spearman rho",
                row_names_gp = gpar(col = label_col),
                column_names_gp = gpar(col = label_col))

pdf("../../out/image/custom/05_heatmap_correlation_GOI_module.pdf",width = 17,height = 12)
draw(ht01)
dev.off()

# significance testing -----------------------------------------------------
# test each gene in the panel against each GOI (Spearman rho + p-value), then FDR-correct across the panel, separately per GOI
test_genes <- setdiff(rownames(average_GOI_RNA_filter),GOI)

df_pval <- purrr::map_dfr(GOI,function(goi){
  purrr::map_dfr(test_genes,function(gene){
    test <- suppressWarnings(
      cor.test(as.numeric(average_GOI_RNA_filter[gene,]),
               as.numeric(average_GOI_RNA_filter[goi,]),
               method = "spearman")
    )
    tibble(Gene = gene,
           target = goi,
           rho = unname(test$estimate),
           pvalue = test$p.value)
  })
}) %>%
  group_by(target) %>%
  mutate(padj = p.adjust(pvalue,method = "BH")) %>%
  ungroup()

# save the full stats table for reporting
df_pval %>%
  arrange(target,padj) %>%
  write_tsv("../../out/table/custom/05_Aletta_correlation_FTH1_GPNMB_stats.tsv")

# wrangle to wide format (rho/padj per GOI) for the scatter
df_pval_wide <- df_pval %>%
  pivot_wider(id_cols = Gene,names_from = target,values_from = c(rho,padj))

# scatter plot focussing only fo the GOIs
df_corr_scatter <- df_pval_wide %>%
  dplyr::rename(FTH1 = rho_FTH1,GPNMB = rho_GPNMB) %>%
  left_join(df_GOI,by = "Gene") %>%
  mutate(sig = case_when(
    padj_FTH1 < 0.05 & padj_GPNMB < 0.05 ~ "both",
    padj_FTH1 < 0.05 ~ "FTH1 only",
    padj_GPNMB < 0.05 ~ "GPNMB only",
    TRUE ~ "ns"
  )) %>%
  mutate(Module = factor(Module))

p_scatter_sig <- df_corr_scatter %>%
  filter(sig != "ns") %>%
  ggplot(aes(x=GPNMB,y=FTH1)) +
  geom_point(aes(col = Module,shape = sig)) +
  ggrepel::geom_text_repel(aes(label = Gene),bg.r = 0.15,bg.color = "white",size = 3) +
  theme_minimal() +
  coord_fixed() +
  geom_hline(yintercept = 0,linetype = "dashed")+
  geom_vline(xintercept = 0,linetype = "dashed") +
  scale_color_manual(values = module_palette,drop = FALSE) +
  labs(title = "Significant genes only")

p_scatter_all <- df_corr_scatter %>%
  ggplot(aes(x=GPNMB,y=FTH1)) +
  geom_point(aes(col = Module)) +
  theme_minimal() +
  coord_fixed() +
  geom_hline(yintercept = 0,linetype = "dashed")+
  geom_vline(xintercept = 0,linetype = "dashed") +
  scale_color_manual(values = module_palette) +
  labs(title = "All panel genes")

p03 <- (p_scatter_all + p_scatter_sig) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "FTH1 vs GPNMB correlation across the gene panel",
                   subtitle = "Spearman rho per gene (pseudobulk per sample, MG cells, all sources pooled)")

ggsave(plot = p03 ,filename = "../../out/image/custom/05_scatter_FTH1_GPNMB_correlation.pdf",width = 15,height = 6)
