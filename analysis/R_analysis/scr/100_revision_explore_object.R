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

# read in the final object ------------------------------------------------
# read in the original dataset
WM_ansinta_nature <- readRDS("../../data/WM_data/WM_data_Absinta_Nature2021/all20_integrated_clean_metadata.rds")
# load the current dataset
scobj <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegAllSoupX_expertAnno.rds")

DimPlot(scobj,group.by = "expertAnno.l1",split.by = "origin",raster=T,label = T)
DimPlot(scobj,group.by = "seurat_clusters",split.by = "origin_alt",raster=T,label = T)
DimPlot(scobj,group.by = "seurat_clusters",split.by = "origin_alt",raster=T,label = T)
DimPlot(scobj,group.by = "seurat_clusters_ref",split.by = "origin_alt",raster=T,label = T)

# check quality of the clusters
df_QC <- scobj@meta.data %>%
  rownames_to_column() %>%
  dplyr::select(rowname,seurat_clusters,nCount_RNA,nFeature_RNA,percent.mt,percent.ribo) %>%
  pivot_longer(names_to = "var",values_to = "value",-c(rowname,seurat_clusters))

df_QC_summary <- df_QC %>%
  group_by(var) %>%
  summarise(med = median(value))

df_QC %>%
  ggplot(aes(x=seurat_clusters,y=value)) +
  geom_violin() +
  geom_boxplot(width = 0.1,outlier.shape = NA) +
  # geom_point(position = position_jitter(width = 0.2),alpha=0.1,size =0.1) +
  facet_wrap(~var,scales="free")+scale_y_sqrt()+
  geom_hline(data = df_QC_summary,aes(yintercept = med),linetype = "dashed",col="red") +
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../../out/image/revision/100_QC_seurat_clusters.pdf",width = 16,height = 8)

# compare with old annotation from martina --------------------------------
# subset only the original data from Martina
test <- subset(scobj,subset = origin_alt == "wm")
test$seurat_clusters_ref <- test$seurat_clusters_ref %>%
  fct_relevel(as.character(0:17))

# plot in comparison with the new embedding
(DimPlot(WM_ansinta_nature,group.by = "seurat_clusters",raster=T,label = T)+ggtitle("original UMAP"))+
(DimPlot(test,group.by = "seurat_clusters_ref",split.by = "origin_alt",raster=T,label = T) + ggtitle("new UMAP"))
ggsave("../../out/image/revision/100_UMAP_wmAbsinta2021_new_vs_old_embedding.pdf",width = 14,height = 6)




