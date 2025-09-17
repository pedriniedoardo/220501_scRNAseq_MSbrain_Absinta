# AIM ---------------------------------------------------------------------
# run pseudobulk analysis following the workflow recommended in the seurat vignette.
# the reference of the test is presented here: https://satijalab.org/seurat/articles/de_vignette.html

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(tidyverse)
library(harmony)
library(ggExtra)
library(ComplexUpset)
library(cowplot)
library(UpSetR)
library(ComplexHeatmap)

# read in the dataset -----------------------------------------------------
# read in the full dataset
sobj <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
DimPlot(sobj,raster = T,group.by = "expertAnno.l1")
DimPlot(sobj,raster = T,group.by = "expertAnno.l1",split.by = "pathology_class")

# test DGE at single cell -------------------------------------------------
# add in one covariate the cell anntation and the stimulation status
sobj$expertAnno.l1.pathology_class <- paste(sobj$expertAnno.l1, sobj$pathology_class, sep = "|")
sobj$expertAnno.l1.pathology_class.sample <- paste(sobj$expertAnno.l1, sobj$pathology_class,sobj$orig.ident, sep = "|")

# set the ident
Idents(sobj) <- "expertAnno.l1.pathology_class"

# run the DGE over the same cell type for stim vs ctrl. The units are the single cells.
# log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
# mono.de <- FindMarkers(sobj, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", verbose = FALSE)

# head(mono.de)
# The p-values obtained from this analysis should be interpreted with caution, because these tests treat each cell as an independent replicate and ignore inherent correlations between cells originating from the same sample. Such analyses have been shown to find a large number of false positive associations, as has been demonstrated by Squair et al., 2021, Zimmerman et al., 2021, Junttila et al., 2022, and others. Below, we show how pseudobulking can be used to account for such within-sample correlation.


# test DGE at pseudobulk --------------------------------------------------
# subset only the cells of interest, in this case only NEU from the CX
sobj_filter <- subset(sobj,subset = origin == "cortex" & (expertAnno.l1 =="EXC NEU"|expertAnno.l1 =="INH NEU"))
table(sobj_filter$orig.ident,sobj_filter$pathology_class,sobj_filter$expertAnno.l1)

# pseudobulk the counts based on donor-condition-celltype
pseudo_sobj <- AggregateExpression(sobj_filter,
                                   assays = "RNA",
                                   slot = "count",
                                   return.seurat = T,
                                   group.by = c("orig.ident","expertAnno.l1","pathology_class"))

# add the covariate for the stimulation per cell type
df_meta_fix <- pseudo_sobj@meta.data %>%
  rownames_to_column() %>%
  mutate(expertAnno.l1 = str_extract_all(rowname,"EXC NEU|INH NEU") %>% unlist(),
         pathology_class = str_extract_all(rowname,"CX_Ctrl|CX_Demye|CX_Mye")%>% unlist()) %>%
  mutate(celltype.stim = paste(expertAnno.l1,pathology_class,sep = "|")) %>%
  column_to_rownames()

pseudo_sobj@meta.data <- df_meta_fix

# explore the dimensionality of the new dataset
pseudo_sobj@meta.data %>%
  group_by(celltype.stim) %>%
  summarise(n = n())

# -------------------------------------------------------------------------
# sample draft to run the pbulk assessment over all the cell types

# pull all the individual cell types
cell_id <- pseudo_sobj@meta.data$expertAnno.l1 %>% unique()

# define the ident for the test
Idents(pseudo_sobj) <- "celltype.stim"

# loop the test over all the cell_ids fro the STIM vs CTRL comparison (per cell id)
# comparisoin 1 mye vs ctrl
# id <- "EXC NEU"
bulk.de.mye <- lapply(cell_id,function(id){
  # keep track the processing
  print(id)
  
  # define the ident
  ident1 <- paste0(id,"|CX_Mye")
  ident2 <- paste0(id,"|CX_Ctrl")
  
  # run the test
  bulk.mono.de <- FindMarkers(object = pseudo_sobj, 
                              ident.1 = ident1, 
                              ident.2 = ident2,
                              test.use = "DESeq2")
  
  # add the gene and cell_id to the data.frame
  bulk.mono.de <- bulk.mono.de %>%
    rownames_to_column("gene") %>%
    mutate(cell_id = id) %>%
    mutate(test = "CX_Mye_vs_CX_Ctrl")
  
  # return the data.frame
  return(bulk.mono.de)
}) %>%
  # make the list as a data.frame
  bind_rows()

# comparisoin 2 demye vs ctrl
bulk.de.demye <- lapply(cell_id,function(id){
  # keep tranck fo the processing
  print(id)
  
  # define the ident
  ident1 <- paste0(id,"|CX_Demye")
  ident2 <- paste0(id,"|CX_Ctrl")
  
  # run the test
  bulk.mono.de <- FindMarkers(object = pseudo_sobj, 
                              ident.1 = ident1, 
                              ident.2 = ident2,
                              test.use = "DESeq2")
  
  # add the gene and cell_id to the data.frame
  bulk.mono.de <- bulk.mono.de %>%
    rownames_to_column("gene") %>%
    mutate(cell_id = id)%>%
    mutate(test = "CX_Demye_vs_CX_Ctrl")
  
  # return the data.frame
  return(bulk.mono.de)
}) %>%
  # make the list as a data.frame
  bind_rows()

# comparisoin 3 demye vs mye
bulk.de.demyeVSmye <- lapply(cell_id,function(id){
  # keep tranck fo the processing
  print(id)
  
  # define the ident
  ident1 <- paste0(id,"|CX_Demye")
  ident2 <- paste0(id,"|CX_Mye")
  
  # run the test
  bulk.mono.de <- FindMarkers(object = pseudo_sobj, 
                              ident.1 = ident1, 
                              ident.2 = ident2,
                              test.use = "DESeq2")
  
  # add the gene and cell_id to the data.frame
  bulk.mono.de <- bulk.mono.de %>%
    rownames_to_column("gene") %>%
    mutate(cell_id = id) %>%
    mutate(test = "CX_Demye_vs_CX_Mye")
  
  # return the data.frame
  return(bulk.mono.de)
}) %>%
  # make the list as a data.frame
  bind_rows()

# merge the three tests
bulk.de <- bind_rows(bulk.de.mye,bulk.de.demye,bulk.de.demyeVSmye)

# save the table
bulk.de %>%
  select(-c(pct.1,pct.2)) %>%
  write_tsv("../../out/table/139_bulk.de.tsv")

# sample plot
bulk.de %>%
  mutate(test_sig = case_when(abs(avg_log2FC)>1 & p_val_adj < 0.05 ~ "sig",
                          T~"not sig")) %>%
  ggplot(aes(x=avg_log2FC,y=-log(p_val_adj))) +
  geom_point(shape = 1,alpha=0.4,aes(col = test_sig)) +
  facet_grid(test~cell_id) +
  geom_vline(xintercept = c(-1,1),col="red",linetype = "dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype = "dashed")+
  theme_bw() +
  scale_color_manual(values = c("black","red")) +
  theme(strip.background = element_blank())

# heatmaps ----------------------------------------------------------------
# pull the scaled values

# design <- model.matrix(~ treat)
# colnames(design)[1] <- c("intercept")
# 
# # save the disign
# saveRDS(design,"../out/object/101_design_mono_pBulk.rds")
# 
# # Create DESeq2 object
# dds <- DESeqDataSetFromMatrix(countData = pseudo_sobj@assays$RNA@counts,
#                               colData = pseudo_sobj@meta.data,
#                               design = design)
# 
# vds_filter <- vst(pseudo_sobj@assays$RNA@counts, blind = F)

mat_filter <- pseudo_sobj@assays$RNA@data %>%
  as.data.frame() %>%
  dplyr::select(contains(c("EXC NEU"))) %>%
  dplyr::select(contains(c("CX_Ctrl","CX_Demye"))) %>%
  as.matrix()

# the DEGs plot
DEG_2 <- bulk.de %>%
  filter(cell_id == "EXC NEU") %>%
  filter(test == "CX_Demye_vs_CX_Ctrl") %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((p_val_adj<0.05)&abs(avg_log2FC)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1) %>%
  pull(gene) %>%
  unique()

# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()

# mat_shr <- mat_filter[rownames(vds_filter) %in% DEG_2, ]
# mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr,useNames = TRUE)
mat_shr <- mat_filter[rownames(mat_filter) %in% DEG_2, ]
mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr,useNames = TRUE)

#
meta_sample <- data.frame(colname = colnames(mat2_shr)) %>%
  left_join(pseudo_sobj@meta.data %>% rownames_to_column("colname"),by=c("colname"))

# make the column of the matrix more readable
# colnames(mat2_shr) <- meta_sample$colname

column_ha_shr <- HeatmapAnnotation(cellid = meta_sample$expertAnno.l1,
                                   treat = meta_sample$pathology_class,
                                   # gender = meta_sample_ENDO$sex,
                                   # facility = meta_sample_ENDO$facility,
                                   col = list(cellid = c("EXC NEU" = "blue",
                                                         "INH NEU" = "yellow"),
                                              treat = c("CX_Ctrl" = "gray",
                                                        "CX_Mye" = "gray30",
                                                        "CX_Demye" = "black")))

ht2_shr <- Heatmap(mat2_shr, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                   name = "exp",
                   column_title = "test data",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr,show_row_names = T,
                   # cluster_rows = F,
                   # right_annotation = row_ha,
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
# pdf("../out/plot/101_heatmap_DEG_plot_pseudobulk_filterExp_shr.pdf",width = 7,height = 14)
draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
# dev.off()


# plot as boxplot
mat2_shr %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene,names_to = "colname",values_to = "scaled_exp") %>%
  left_join(meta_sample,by = "colname") %>%
  ggplot(aes(x=pathology_class,y=scaled_exp)) + geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(width = 0.1))+facet_wrap(expertAnno.l1~gene)+theme_bw()+theme(strip.background = element_blank())


# -------------------------------------------------------------------------
# explore some genes that are labelled as DEGs in both versions of the analysis
# top 5
top5_common <- 
  bulk.de %>%
  filter(cell_id == "EXC NEU") %>%
  filter(test == "CX_Demye_vs_CX_Ctrl") %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((p_val_adj<0.05)&abs(avg_log2FC)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1) %>%
  dplyr::slice(1:5) %>%
  pull(gene)

# plot top5 by condition
Idents(sobj) <- "expertAnno.l1.pathology_class"
VlnPlot(sobj, features = top5_common, idents = c("EXC NEU|CX_Ctrl", "EXC NEU|CX_Demye"),
        group.by = "expertAnno.l1.pathology_class.sample",ncol = 1,
        raster = F)
