# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SeuratData)
library(ggridges)
library(ComplexHeatmap)
library(cowplot)

# read in the data --------------------------------------------------------
# in this case I want to use the Leng 2021 dataset as reference to map the cells of martina's dataset
# read in the reference dataset for the leng dataset
ref_WM_original <- readRDS("../../data/all20_immune.rds")
ref_WM <- readRDS("../../data/all20_immune_model.rds")

# add the original cluster id to the new immune subset with the model
LUT_clusterID <- left_join(
  ref_WM@meta.data %>%
    rownames_to_column("barcode") %>%
    dplyr::select(barcode,seurat_clusters_new = seurat_clusters),
  ref_WM_original@meta.data %>%
    rownames_to_column("barcode") %>%
    dplyr::select(barcode,seurat_clusters_original = seurat_clusters),by = "barcode") %>%
  mutate(seurat_clusters_new = paste0("new_",seurat_clusters_new),
         seurat_clusters_original = paste0("org_",seurat_clusters_original))

table(LUT_clusterID$seurat_clusters_new,
      LUT_clusterID$seurat_clusters_original)

DimPlot(ref_WM,label = T)
DimPlot(ref_WM_original,label = T)

# update the metadata
ref_WM$seurat_clusters_new <- LUT_clusterID$seurat_clusters_new
ref_WM$seurat_clusters_original <- LUT_clusterID$seurat_clusters_original

#
DimPlot(ref_WM,group.by = "seurat_clusters_original",label=T)
# # get the metadata from the other object
# meta <- ref_WM@meta.data %>%
#   rownames_to_column(var = "barcodes") %>%
#   # add the macroclassification
#   mutate(clusterCellType = case_when(seurat_clusters %in% c(0,1,2,3,6) ~"OLIGO",
#                                      seurat_clusters %in% c(7,15) ~"NEU",
#                                      seurat_clusters %in% c(8) ~"OPC",
#                                      seurat_clusters %in% c(11,13) ~"VAS",
#                                      seurat_clusters %in% c(16) ~"LYM",
#                                      seurat_clusters %in% c(5,10,17) ~"IMM",
#                                      seurat_clusters %in% c(4,9,12,14) ~"AST"))
# 
# ref_WM$clusterCellType <- meta$clusterCellType
# 
# ref_WM@meta.data
# print the top markers per cluster
DefaultAssay(ref_WM)<-"RNA"
Idents(ref_WM) <- "seurat_clusters_original"

# data.combined.markers <- FindAllMarkers(ref_WM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# # save the table of all markers
# data.combined.markers %>%
#   write_tsv("out/table/FindAllMarkers_WMimm.tsv")
# 
# # top 100 per cluster
# data.combined.markers %>%
#   group_by(cluster) %>%
#   mutate(rank = rank(order(p_val_adj, -abs(avg_log2FC)), ties.method='first')) %>%
#   arrange(cluster,rank) %>%
#   filter(rank < 101) %>%
#   write_tsv("out/table/FindAllMarkers_WMimm_top100.tsv")

# read in the query dataset
query <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_MG_harmonySkipIntegration.rds")

# fix the meta to add the general cell annotation
query$seurat_clusters_fix1 <- paste0("clu_",query$seurat_clusters)

query@meta.data

DimPlot(query,group.by = "seurat_clusters_fix1",label = T)
# -------------------------------------------------------------------------
# to work the reference should be the integrated dataset
# DefaultAssay(ref_WM)<-"integrated"
DefaultAssay(ref_WM)<-"integrated"
# the RNA slot do not work due to the missing variable features
# DefaultAssay(ref_WM)<-"RNA"
dim(ref_WM)
ref_WM@meta.data

# in the vignetted the query isn the RNA slot
# DefaultAssay(query) <- "integrated"
DefaultAssay(query) <- "RNA"
dim(query)

test.anchors_0 <- FindTransferAnchors(reference = ref_WM,
                                      query = query,
                                      dims = 1:30,
                                      reference.reduction = "pca")

predictions_0 <- TransferData(anchorset = test.anchors_0,
                              refdata = ref_WM$seurat_clusters_original,
                              dims = 1:30)

# show the distribution of the max scores
predictions_0 %>%
  ggplot(aes(x=prediction.score.max))+geom_histogram()+theme_bw()
ggsave("../../out/image/histo_prediction_label_transfer_refWMImmune_queryWMCX.pdf",width = 5,height = 4)

dim(predictions_0)
# add the predictions ot the query dataset
query_transfer <- AddMetaData(query, metadata = predictions_0)

# -------------------------------------------------------------------------
# plot on the original query dataset
# DimPlot(query, reduction = "umap", group.by = "predicted.id")
# see how is different from the integrated one
# DimPlot(ref_WM, reduction = "umap", group.by = "clusterCellType")

# add the meta to the coordinates
data_transfer <- left_join(
  # get the coordinates of the query dataset
  query_transfer@reductions$umap@cell.embeddings %>%
    data.frame() %>%
    rownames_to_column(var = "barcodes"),
  
  # get the metadata from the query dataset, after adding the prediciton information. use a threshold for the score
  query_transfer@meta.data %>%
    rownames_to_column(var = "barcodes") %>%
    # this is based on martina's table, I need to ask her how she provided the accuracy of the imputation
    mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                    T ~ "uncertain"))
  ,"barcodes")

# notice that basically all the oligos in the braak0 dataset are not correctly mapped 
table(data_transfer$robust_score)

# -------------------------------------------------------------------------
# average the position of the clusters
data_transfer_avg_raw <- data_transfer %>% group_by(predicted.id) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# plot the raw prediction based on the reference dataset
ggplot(label= TRUE)+
  geom_point(data = data_transfer,aes(x = UMAP_1,y = UMAP_2, col = predicted.id),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_transfer_avg_raw,aes(x = UMAP_1,y = UMAP_2,label = predicted.id),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("../../out/image/UMAP_label_transfer_refWMImmune_queryWMCX_raw.pdf",width = 5,height = 4)

# plot the raw prediction based on the reference dataset robust
data_transfer_unc <- data_transfer %>% 
  mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                  T ~ "uncertain")) %>% 
  filter(robust_score == "uncertain")
#
data_transfer_defined <- data_transfer %>% 
  mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                  T ~ "uncertain")) %>% 
  filter(robust_score != "uncertain")

data_transfer_avg_robust <- data_transfer_defined %>% group_by(predicted.id) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

ggplot(label= TRUE)+
  geom_point(data = data_transfer_unc,aes(x = UMAP_1,y = UMAP_2),size=0.1,alpha=0.1,col="gray") +
  geom_point(data = data_transfer_defined,aes(x = UMAP_1,y = UMAP_2, col = predicted.id),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_transfer_avg_robust,aes(x = UMAP_1,y = UMAP_2,label = predicted.id),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("../../out/image/UMAP_label_transfer_refWMImmune_queryWMCX_robust.pdf",width = 5,height = 4)

# -------------------------------------------------------------------------
# how many cells of the quey astro are high confidence astro in the reference dataset
# query_transfer@meta.data %>%
#   # pull the astro from the query annotation
#   filter(seurat_clusters_fix1 %in% c("AST")) %>%
#   mutate(robust_score = case_when(prediction.score.max>0.7~predicted.id,
#                                   T ~ "uncertain")) %>%
#   summarise(tot = n(),
#             prop_correct = mean(robust_score != "uncertain"))
# #
# query_braak6@meta.data %>%
#   # pull the astro from the query annotation
#   filter(seurat_clusters_fix1 %in% c("AST")) %>%
#   mutate(robust_score = case_when(prediction.score.max>0.7~predicted.id,
#                                   T ~ "uncertain")) %>%
#   summarise(tot = n(),
#             prop_correct = mean(robust_score != "uncertain"))

# plot the predicted score per per cell type
# since I am placing my cells on a reference it makes sense to me to show the score of how well each of my cells i mappe on the reference
id_factor <- query_transfer@meta.data %>% 
  group_by(seurat_clusters_fix1) %>% 
  summarise(med = median(prediction.score.max)) %>% 
  arrange(med) %>% 
  pull(seurat_clusters_fix1)

# plot the mapping scores for the original label ID
query_transfer@meta.data %>%
  mutate(seurat_clusters_fix1= factor(seurat_clusters_fix1,levels = id_factor)) %>% 
  ggplot(aes(x=prediction.score.max,y=seurat_clusters_fix1)) +
  geom_density_ridges(scale = 0.9) +
  theme_bw() +
  geom_vline(xintercept = 0.70,linetype="dashed",col="red",alpha=0.5)
ggsave("../../out/image/histo_label_transfer_refWMImmune_queryWMCX.pdf",width = 3,height = 4)

# plot the mapping scores for the ref labels
id_factor2 <- query_transfer@meta.data %>% 
  group_by(predicted.id) %>% 
  summarise(med = median(prediction.score.max)) %>% 
  arrange(med) %>% 
  pull(predicted.id)

query_transfer@meta.data %>%
  mutate(predicted.id= factor(predicted.id,levels = id_factor2)) %>% 
  ggplot(aes(x=prediction.score.max,y=predicted.id)) +
  geom_density_ridges(scale = 0.9) +
  theme_bw() +
  geom_vline(xintercept = 0.70,linetype="dashed",col="red",alpha=0.5)
ggsave("../../out/image/histo_label_transfer_refWMImmune_queryWMCX2.pdf",width = 3,height = 4)

# heatmap for the assignament ---------------------------------------------
prop_table_tot <- query_transfer@meta.data %>%
  mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                  T ~ "uncertain")) %>% 
  # filter(robust_score_subclass != "uncertain") %>%
  # filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters_fix1,predicted.id) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters_fix1,predicted.id,prop) %>%
  pivot_wider(names_from = seurat_clusters_fix1,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.id") %>%
  as.matrix()

prop_table_robust <- query_transfer@meta.data %>%
  mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                  T ~ "uncertain")) %>% 
  filter(robust_score != "uncertain") %>%
  group_by(seurat_clusters_fix1,predicted.id) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters_fix1,predicted.id,prop) %>%
  pivot_wider(names_from = seurat_clusters_fix1,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.id") %>%
  as.matrix()

pdf("../../out/image/heatmap_label_transfer_refWMImmune_queryWMCX_raw.pdf",height = 3,width = 3)
Heatmap(prop_table_tot,
        name = "prop", 
        column_title = "predicted WM id raw",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

pdf("../../out/image/heatmap_label_transfer_refWMImmune_queryWMCX_robust.pdf",height = 3,width = 3)
Heatmap(prop_table_robust,
        name = "prop", 
        column_title = "predicted WM id robust",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# Unimodal UMAP Projection ------------------------------------------------
# In Seurat v4, we also enable projection of a query onto the reference UMAP structure. This can be achieved by computing the reference UMAP model and then calling MapQuery() instead of TransferData().
# this one run on the subset will overwrite the model
# if I don't run it the result fo the mapping looks weird
# ref_WM <- RunUMAP(ref_WM,
#                   dims = 1:30,
#                   reduction = "pca",
#                   return.model = TRUE)
query_transfer <- MapQuery(anchorset = test.anchors_0,
                           reference = ref_WM,
                           query = query_transfer,
                           refdata = list(celltype = "seurat_clusters_original"),
                           reference.reduction = "pca", reduction.model = "umap")
#
# We can now visualize the query cells alongside our reference.
# -------------------------------------------------------------------------
#
# DimPlot(ref_WM, reduction = "umap", group.by = "clusterCellType", label = TRUE, label.size = 3,repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
# DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

# add the meta to the coordinates
data2_transfer <- left_join(
  # get the coordinates of the query dataset
  query_transfer@reductions$ref.umap@cell.embeddings %>%
    data.frame() %>%
    rownames_to_column(var = "barcodes"),
  
  # get the metadata from the query dataset, after adding the prediciton information. use a threshold for the score
  query_transfer@meta.data %>%
    rownames_to_column(var = "barcodes") %>%
    # this is based on martina's table, I need to ask her how she provided the accuracy of the imputation
    mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                    T ~ "uncertain"))
  ,"barcodes")

# average the position of the clusters
data2_transfer_avg_raw <- data2_transfer %>% group_by(seurat_clusters_fix1) %>% select(refUMAP_1, refUMAP_2) %>% summarize_all(mean)

# plot the raw prediction based on the reference dataset
# save the reference UMAPS
UMAP_ref_WM <- DimPlot(ref_WM,group.by = "seurat_clusters_original",label = T)
UMAP_ref_WM$data %>%
  ggplot(label= TRUE)+
  UMAP_ref_WM$layers[[1]]+
  UMAP_ref_WM$layers[[2]]+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("../../out/image/UMAP_ref_WM.pdf",width = 6,height = 4)

# save the image with no label and no legend
UMAP_ref_WM$data %>%
  ggplot(label= F)+
  UMAP_ref_WM$layers[[1]]+
  # UMAP_ref_WM$layers[[2]]+
  theme_bw()+
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("../../out/image/UMAP_ref_WM_no_lab.pdf",width = 4,height = 4)

ggplot(label= TRUE)+
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  geom_point(data = data2_transfer,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data2_transfer_avg_raw,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))
ggsave("../../out/image/UMAP_label_transfer_refWMImmune_queryWMCX_raw_inverse.pdf",width = 6,height = 4)

# save the image with no label and no legend
ggplot(label= F)+
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  geom_point(data = data2_transfer,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = data2_transfer_avg_raw,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+
  theme(legend.position = "none")
ggsave("../../out/image/UMAP_label_transfer_refWMImmune_queryWMCX_raw_inverse_no_lab.pdf",width = 5,height = 4)

data2_transfer %>% 
  filter(robust_score != "uncertain") %>% 
  group_by(seurat_clusters_fix1,robust_score) %>% 
  summarise(n = n()) %>% 
  print(n=70)

# divide the dataset into uncertain and not
data2_transfer_unc <- data2_transfer %>%
  filter(robust_score == "uncertain")
#
data2_transfer_defined <- data2_transfer %>%
  filter(robust_score != "uncertain")
# # focus only on cluster with more than 20 cells with robust prediction
# group_by(seurat_clusters_fix1,robust_score) %>% 
# mutate(tot_n = n()) %>% 
# filter(tot_n>20)

# # focus only on cluster with more than 20 cells with robust prediction
# data2_transfer_defined  %>% 
#   group_by(seurat_clusters_fix1,robust_score) %>% 
#   summarise(tot_n = n()) %>% 
#   filter(seurat_clusters=="cluster_10")

# average the position of the clusters
data2_transfer_avg <- data2_transfer_defined %>% group_by(seurat_clusters_fix1) %>% select(refUMAP_1, refUMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE) +
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  # geom_point(data = data2_transfer_unc,aes(x = refUMAP_1,y = refUMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data2_transfer_defined,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data2_transfer_avg,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))
ggsave("../../out/image/UMAP_label_transfer_refWMImmune_queryWMCX_robust_inverse.pdf",width = 6,height = 4)

# print also without the labels
ggplot(label= F) +
  geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  # geom_point(data = data2_transfer_unc,aes(x = refUMAP_1,y = refUMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data2_transfer_defined,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = data2_transfer_avg,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+
  theme(legend.position = "none")
ggsave("../../out/image/UMAP_label_transfer_refWMImmune_queryWMCX_robust_inverse_no_lab.pdf",width = 5,height = 4)

# ggplot(label= TRUE) +
#   geom_point(data = ref_WM@reductions$umap@cell.embeddings %>%
#                data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
#   # geom_point(data = data2_transfer_unc,aes(x = refUMAP_1,y = refUMAP_2),size=0.3,alpha=0.1,col="gray") +
#   geom_point(data = data2_transfer_defined,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
#   # labs(color= "Clusters") +
#   ggrepel::geom_text_repel(data = data2_transfer_avg,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+theme_bw()+
#   guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+ facet_wrap(~treat)+theme(strip.background = element_blank())
# 
# df1 <- data2_transfer %>% 
#   group_by(treat) %>% 
#   summarise(count_tot = n())
# 
# df2 <- data2_transfer_defined %>% 
#   group_by(treat) %>% 
#   summarise(count_robust = n())
# 
# left_join(df1,df2) %>% 
#   mutate(prop = count_robust/count_tot)
# 
