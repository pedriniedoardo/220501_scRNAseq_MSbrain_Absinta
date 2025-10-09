# AIM ---------------------------------------------------------------------
# try simple correaltion analysis for exploring the similarity between out sc dataset and Eliana's dataset
# load the batch-corrected, normalized expression matrix (assay(vst_CellidSmall_VAS_batch)), and calculate the correlation between the pure bulk endothelial samples and each of the pseudo-bulk cell types from the brain dataset
# try not to use only the subset of vas genes

# libraries ---------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(DESeq2)

# read in the data --------------------------------------------------------
# load the scaled, batch corrected matrix from VAS specific genes.
vst_batch <- readRDS("../../out/object/204_vst_filterBatch_Eliana_Absinta_CellidSmall.rds")
# vst_batch_VAS <- readRDS("../../out/object/204_vst_filterBatch_Eliana_Absinta_CellidSmall.rds")

# processing --------------------------------------------------------------
# produce the matrix and metadata
vst_batch_mat <- assay(vst_batch)
metadata <- colData(vst_batch)

# Use the cor() function to compute the correlation (Spearman is often preferred over Pearson as it's rank-based and less sensitive to outliers) between every sample.
correlation_matrix <- cor(vst_batch_mat,method = "spearman")

# Visualize the results as a heatmap
pdf("../../out/image/205_heatmap_all_correlation_batch.pdf",width = 7,height = 6)
Heatmap(correlation_matrix,col = viridis::viridis(option = "turbo",n = 10))
dev.off()

# Use the cor() function to compute the correlation (Spearman is often preferred over Pearson as it's rank-based and less sensitive to outliers) between every bulk sample and every pseudo-bulk profile.
bulk_samples <- rownames(metadata[metadata$batch == "Eliana", ])
sc_samples <- rownames(metadata[metadata$batch == "Absinta", ])

# Calculate Spearman correlation
correlation_matrix_VAS2 <- cor(vst_batch_mat[, bulk_samples], 
                               vst_batch_mat[, sc_samples], 
                               method = "spearman")

# Visualize the results as a heatmap
pdf("../../out/image/205_heatmap_short_correlation_batch.pdf",width = 7,height = 6)
Heatmap(correlation_matrix_VAS2, col = viridis::viridis(option = "turbo",n = 10))
dev.off()
