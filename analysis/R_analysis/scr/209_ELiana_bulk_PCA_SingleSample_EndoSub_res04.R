# AIM ---------------------------------------------------------------------
# put together the bulk RNAseq data from Eliana and out pBulk data from the brain dataset
# in this case combine all the sample from Eliana in just one sample
# in this case I am aggregating the samples from the endothelia subclusters only

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)
library(DESeq2)
library(RNAseqQC)
library(limma)
library(ashr)
library(magick)
library(UpSetR)
library(ComplexHeatmap)
library(pals)
library(scales)
library(patchwork)
library(enrichR)
library(GGally)
library(presto)
library(SeuratWrappers)

# produce the table for the different runs --------------------------------
# read in the data from the sc dataset
scobj <- readRDS("../../out/object/130_VAS_subcluster_HarmonySample.rds")
DimPlot(scobj,raster = F,group.by = "RNA_snn_res.0.4",label = T)

# simple aggregation of the data by the specific resolution
cts_sample_CellidSmall <- AggregateExpression(object = scobj,
                                              group.by = c("RNA_snn_res.0.4"),
                                              assays = 'RNA',
                                              slot = "counts",
                                              return.seurat = FALSE)

counts_CellidSmall <- cts_sample_CellidSmall$RNA
dim(counts_CellidSmall)

# to avoid issues, make the column name for the subclusters, a bit better
colnames(counts_CellidSmall) <- paste0("clu_",colnames(counts_CellidSmall))

# remove the sc object not needed anymore
# remove(scobj)
# gc()

# read in the expresson data from Eliana bulk rnaseq data from ECFC cells
# dds <- readRDS(file = "../../out/object/201_dds_all.rds")
dds <- readRDS(file = "../../out/object/201_dds_all_update.rds")

count_ecfc_all <- counts(dds)
meta_ecfc_all <- colData(dds) %>%
  data.frame()

# aggregate the reads and the metadata in one sample only
count_ecfc <- rowSums(count_ecfc_all) %>%
  data.frame(ECFC = .)

meta_ecfc <- data.frame(sample = "ECFC",Group = "ECFC") %>%
  mutate(rowname = sample) %>%
  column_to_rownames()

head(count_ecfc)
meta_ecfc

# wrangling ---------------------------------------------------------------
# join the two tables make an inner join to so drop the non common genes across the two datasets
counts_full_CellidSmall <- inner_join(counts_CellidSmall %>%
                                        data.frame() %>%
                                        rownames_to_column("gene"),
                                      count_ecfc %>%
                                        data.frame() %>%
                                        rownames_to_column("gene"),by = "gene") %>%
  column_to_rownames("gene")

# compare the dimensions
dim(counts_full_CellidSmall)
dim(count_ecfc)
dim(counts_CellidSmall)

# 2. generate sample level metadata
# metadata from the aggregated sc dataset
meta_CellidSmall_short <- data.frame(sample = colnames(counts_CellidSmall) %>% unname()) %>%
  mutate(Group = "mix") %>%
  mutate(batch = "Absinta") %>%
  mutate(rowname = str_replace_all(sample,pattern = "\\s+","\\.")) %>%
  mutate(cellid = sample)

# metadata fro the bulk ecfc dataset
meta_ecfc_short <- meta_ecfc %>%
  dplyr::select(sample,Group) %>%
  mutate(batch = "Eliana") %>%
  remove_rownames() %>%
  mutate(rowname = str_replace_all(sample,pattern = "-","\\.")) %>%
  mutate(cellid = "ECFC") %>%
  mutate(test = "ECFC")

# build the final metadata
colData_full_CellidSmall <- bind_rows(meta_CellidSmall_short,meta_ecfc_short) %>%
  column_to_rownames("rowname") %>%
  .[colnames(counts_full_CellidSmall),]

# save the metadata and the full matrix
saveRDS(counts_full_CellidSmall,"../../out/object/209_count_Eliana_Absinta_CellidSmall_EndoSub_res04.rds")
saveRDS(colData_full_CellidSmall,"../../out/object/209_meta_Eliana_Absinta_CellidSmall_EndoSub_res04.rds")

# bulkd the DESeq2 object -------------------------------------------------
# build the model, I am not really using it.
batch <- colData_full_CellidSmall$batch
design_CellidSmall <- model.matrix(~ batch)
colnames(design_CellidSmall)[1] <- c("intercept")

# save the disign
saveRDS(design_CellidSmall,"../../out/object/209_design_ELiana_Absinta_CellidSmall_EndoSub_res04.rds")

# Create DESeq2 object   
dds_CellidSmall <- DESeqDataSetFromMatrix(countData = counts_full_CellidSmall,
                                          colData = colData_full_CellidSmall,
                                          design = design_CellidSmall)

# filter
# filter out lowly expressed genes
keep_features <- edgeR::filterByExpr(counts(dds_CellidSmall), group = colData_full_CellidSmall$batch)
dds_CellidSmall_filter <- dds_CellidSmall[keep_features,]

# check the dimenstion before and after filtering
dim(dds_CellidSmall)
dim(dds_CellidSmall_filter)

# scale the data using vst. blind = T
vst_CellidSmall_filter <- vst(dds_CellidSmall_filter, blind = T)

# save the objecs
saveRDS(dds_CellidSmall_filter,"../../out/object/209_dds_filter_Eliana_Absinta_CellidSmall_EndoSub_res04.rds")
saveRDS(vst_CellidSmall_filter,"../../out/object/209_vst_filter_Eliana_Absinta_CellidSmall_EndoSub_res04.rds")

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/image/209_heatmap_Eliana_Absinta_CellidSmall_EndoSub_res04.pdf",width = 12,height = 8)
set.seed(26)
hm <- plot_sample_clustering(vst_CellidSmall_filter,
                             anno_vars = c("batch","Group","cellid"),
                             distance = "euclidean")
draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 40), "mm"))
dev.off()

# pull the actual distances
hm@matrix

# # confirm the calculation of the distance
# # this is using all the features
# mat <- assay(vst_CellidSmall_filter)
# head(mat)
# dist(t(mat))
# 
# # uisng the top 500 most variable features
# # pull the variance per gene
# matrixStats::rowVars(mat) %>% head()
# 
# # define the matrix
# mat_top500 <- mat[matrixStats::rowVars(mat) %>%
#                     order(decreasing = T) %>%
#                     head(500), ]
# 
# head(mat_top500)
# str(mat_top500)
# 
# dist(t(mat_top500))

# pull only the ECFC comparison
mat_dist_ECFC <- hm@matrix %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  dplyr::select(ID,ECFC) %>%
  filter(ID != "ECFC") %>%
  column_to_rownames("ID") 


ht <- Heatmap(mat_dist_ECFC,
              show_column_names = T,
              raster_by_magick = T
              ,show_row_dend = T, use_raster = T,
              name = "eucledian \n distance",
              # column_title = "sample",
              col = viridis::viridis(option = "plasma",n = 10),
              
              # row_names_gp = gpar(fontsize = 3),
              # top_annotation = column_meta_sample_prop,
              show_row_names = T
              # cluster_rows = F,
              # right_annotation = row_ha,
              # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
pdf("../../out/image/209_heatmapDistanceECFC_Eliana_Absinta_CellidSmall_EndoSub_res04.pdf",width = 5,height = 4)
draw(ht,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2,2,2,60), "mm"))
dev.off()

# PCA ---------------------------------------------------------------------
# make pca over the non batch corrected values
plot_vst_CellidSmall <- plotPCA(vst_CellidSmall_filter,
                                intgroup = c("batch","Group","cellid","test")) +
  theme_bw()

# extract the pca values from batch uncorrected value of expression
df_plot_PCA <- plot_vst_CellidSmall$data

# fix the colors, set the costum colors palette for the cell ids
myColors <- alphabet(length(unique(df_plot_PCA$cellid)))
names(myColors) <- unique(df_plot_PCA$cellid)
show_col(myColors)

custom_scale_color <- scale_color_manual(values = myColors)

# plot PCA using the new palette of colors
# all sample no batch correction
p1 <- df_plot_PCA %>%
  # make CSF.MS a single category
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=cellid)) +
  geom_point(size =4,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vst_CellidSmall$labels[1]) + xlab(plot_vst_CellidSmall$labels[2]) +
  custom_scale_color +
  ggtitle("All")

# make one coloring by batch
p2 <- df_plot_PCA %>%
  # make CSF.MS a single category
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=batch)) +
  geom_point(size =4,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vst_CellidSmall$labels[1]) + xlab(plot_vst_CellidSmall$labels[2])+
  ggtitle("Batch")

# make one coloring by test
p3 <- df_plot_PCA %>%
  # make CSF.MS a single category
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=test)) +
  geom_point(size =4,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vst_CellidSmall$labels[1]) + xlab(plot_vst_CellidSmall$labels[2])+
  ggtitle("Endo")

# plot the panel
p1 + p2 + p3 + plot_annotation("No batch correction")
ggsave("../../out/image/209_PCA_Eliana_Absinta_CellidSmall_EndoSub_res04.pdf",width = 20,height = 6)

# explore more pcs trend --------------------------------------------------
# pull more PC
rv <- rowVars(assay(vst_CellidSmall_filter),useNames = TRUE)
# select the ntop genes by variance
select_var <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
test <- prcomp(t(assay(vst_CellidSmall_filter)[select_var,]))$x %>% 
  data.frame() %>% 
  rownames_to_column("sample")

# plot more PC by condition
left_join(plot_vst_CellidSmall$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample")) %>%
  # ggpairs(columns = 5:14,ggplot2::aes(colour=condition),upper = "blank")+
  ggpairs(columns = 7:15,ggplot2::aes(colour=test),upper = "blank",legend = 1) +
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# scale_fill_manual(values = myColors) +
# scale_color_manual(values = myColors)
ggsave("../../out/image/209_panel_PC_Eliana_Absinta_CellidSmall_EndoSub_res04.pdf",width = 20,height = 20)

# explore pc score by metadata fo the samples
test2 <- left_join(plot_vst_CellidSmall$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample"))

test_df1 <- test2 %>%
  dplyr::select(batch,Group,cellid,test,name) %>%
  # mutate(exposure = factor(exposure,levels = c("0","6","24"),labels = c("00","06","24"))) %>%
  # dplyr::rename(approx.time = Approx.Time.between.collection.and.processing) |> 
  pivot_longer(names_to = "var_1",values_to = "value_1",-c(name))

test_df2 <- test2 %>%
  dplyr::select(name,PC1:PC9) %>%
  pivot_longer(names_to = "var_2",values_to = "value_2",-c(name))

left_join(test_df1,test_df2,by=c("name")) %>%
  mutate(comparison = paste0(var_1,"_vs_",var_2)) %>%
  ggplot(aes(x=value_1,y=value_2)) +
  facet_wrap(~comparison,scales = "free",ncol=9) +
  # facet_grid(var_1~var_2,scales = "free_x",drop = T) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/image/209_panel_metadata_PC_Eliana_Absinta_CellidSmall_EndoSub_res04.pdf",width = 20,height = 9)

# batch removal -----------------------------------------------------------
# try to remove the batch from the scaled expression
vst_CellidSmall_filter_batch <- vst_CellidSmall_filter

mat <- assay(vst_CellidSmall_filter)
mat_batch <- limma::removeBatchEffect(mat, vst_CellidSmall_filter$batch)
assay(vst_CellidSmall_filter_batch) <- mat_batch

# save the object after batch correction
saveRDS(vst_CellidSmall_filter_batch,"../../out/object/209_vst_filterBatch_Eliana_Absinta_CellidSmall_EndoSub_res04.rds")

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/image/209_heatmap_Eliana_Absinta_CellidSmall_batch_EndoSub_res04.pdf",width = 12,height = 8)
set.seed(26)
hm_batch <- plot_sample_clustering(vst_CellidSmall_filter_batch,
                                   anno_vars = c("batch","Group","cellid"),
                                   distance = "euclidean")
draw(hm_batch,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 40), "mm"))
dev.off()

# pull the actual distances
hm_batch@matrix

# pull only the ECFC comparison
mat_dist_ECFC_batch <- hm_batch@matrix %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  dplyr::select(ID,ECFC) %>%
  filter(ID != "ECFC") %>%
  column_to_rownames("ID") 

ht_batch <- Heatmap(mat_dist_ECFC_batch,
              show_column_names = T,
              raster_by_magick = T
              ,show_row_dend = T, use_raster = T,
              name = "eucledian \n distance",
              # column_title = "sample",
              col = viridis::viridis(option = "plasma",n = 10),
              
              # row_names_gp = gpar(fontsize = 3),
              # top_annotation = column_meta_sample_prop,
              show_row_names = T
              # cluster_rows = F,
              # right_annotation = row_ha,
              # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
pdf("../../out/image/209_heatmapDistanceECFC_Eliana_Absinta_CellidSmall_batch_EndoSub_res04.pdf",width = 5,height = 4)
draw(ht_batch,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2,2,2,60), "mm"))
dev.off()

# PCA ---------------------------------------------------------------------
plot_vst_CellidSmall_batch <- plotPCA(vst_CellidSmall_filter_batch,
                                      intgroup = c("batch","Group","cellid","test")) +
  theme_bw()

# extract the pca values
df_plot_PCA_batch <- plot_vst_CellidSmall_batch$data

# fix the colors
# set the costum colors
# myColors <- alphabet(length(unique(df_plot_PCA$cellid)))
# names(myColors) <- unique(df_plot_PCA$cellid)
# show_col(myColors)
# 
# custom_scale_color <- scale_color_manual(values = myColors)

# all sample after batch correction
p1_batch <- df_plot_PCA_batch %>%
  # make CSF.MS a single category
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=cellid)) +
  geom_point(size =4,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vst_CellidSmall_batch$labels[1]) + xlab(plot_vst_CellidSmall_batch$labels[2]) +
  custom_scale_color +
  ggtitle("All")

# make one coloring by batch
p2_batch <- df_plot_PCA_batch %>%
  # make CSF.MS a single category
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=batch)) +
  geom_point(size =4,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vst_CellidSmall_batch$labels[1]) + xlab(plot_vst_CellidSmall_batch$labels[2])+
  ggtitle("Batch")

# make one coloring by test
p3_batch <- df_plot_PCA_batch %>%
  # make CSF.MS a single category
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  ggplot(aes(x=PC1,y=PC2,col=test)) +
  geom_point(size =4,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vst_CellidSmall$labels[1]) + xlab(plot_vst_CellidSmall$labels[2])+
  ggtitle("Endo")


# plot the panel
p1_batch + p2_batch + p3_batch + plot_annotation("batch correction")
ggsave("../../out/image/209_PCA_Eliana_Absinta_CellidSmall_batch_EndoSub_res04.pdf",width = 20,height = 6)

# explore more pcs trend --------------------------------------------------
# pull more PC
rv_batch <- rowVars(assay(vst_CellidSmall_filter_batch),useNames = TRUE)
# select the ntop genes by variance
select_var_batch <- order(rv_batch, decreasing=TRUE)[seq_len(min(500, length(rv_batch)))]
# perform a PCA on the data in assay(x) for the selected genes
test_batch <- prcomp(t(assay(vst_CellidSmall_filter_batch)[select_var_batch,]))$x %>% 
  data.frame() %>% 
  rownames_to_column("sample")

# plot more PC by condition
left_join(plot_vst_CellidSmall_batch$data %>% dplyr::select(-c("PC1","PC2")),test_batch,by = c("name"="sample")) %>%
  # ggpairs(columns = 5:14,ggplot2::aes(colour=condition),upper = "blank")+
  ggpairs(columns = 7:15,ggplot2::aes(colour=test),upper = "blank",legend = 1) +
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# scale_fill_manual(values = myColors) +
# scale_color_manual(values = myColors)
# ggsave("../../out/image/209_panel_PC_Eliana_Absinta_CellidSmall_batch_EndoSub_res04.pdf",width = 20,height = 20)

# explore pc score by metadata fo the samples
test2_batch <- left_join(plot_vst_CellidSmall_batch$data %>% dplyr::select(-c("PC1","PC2")),test_batch,by = c("name"="sample"))

test_df1_batch <- test2_batch %>%
  dplyr::select(batch,Group,cellid,name,test) %>%
  # mutate(exposure = factor(exposure,levels = c("0","6","24"),labels = c("00","06","24"))) %>%
  # dplyr::rename(approx.time = Approx.Time.between.collection.and.processing) |> 
  pivot_longer(names_to = "var_1",values_to = "value_1",-c(name))

test_df2_batch <- test2_batch %>%
  dplyr::select(name,PC1:PC9) %>%
  pivot_longer(names_to = "var_2",values_to = "value_2",-c(name))

left_join(test_df1_batch,test_df2_batch,by=c("name")) %>%
  mutate(comparison = paste0(var_1,"_vs_",var_2)) %>%
  ggplot(aes(x=value_1,y=value_2)) +
  facet_wrap(~comparison,scales = "free",ncol=9) +
  # facet_grid(var_1~var_2,scales = "free_x",drop = T) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/image/209_panel_metadata_PC_Eliana_Absinta_CellidSmall_batch_EndoSub_res04.pdf",width = 20,height = 9)
