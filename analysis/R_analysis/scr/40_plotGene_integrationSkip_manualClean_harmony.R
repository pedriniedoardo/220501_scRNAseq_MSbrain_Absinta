# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(enrichR)
library(patchwork)
library(cowplot)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegAllSoupX_expertAnno.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")

# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("TSPO")

table(data.combined@meta.data$annotation_confident)

# add a broader cell grouping
# data.combined@meta.data$cell_type2 <- data.combined@meta.data |> 
#   mutate(cell_type2 = case_when(# according to label transfer cluster 12 and 6 are Bipolar cells
#     cell_type %in% c("BP_0","BP_10","BP_11","BP_13","BP_16","BP_17","BP_5","BP_8","12","6")~"BP",
#     T~cell_type)) |> 
#   pull(cell_type2)

# generate the table for the plots ----------------------------------------
# get the metadata from the other object
meta <- data.combined@meta.data %>%
  rownames_to_column(var = "barcodes")

# extrac the expression value
df_exp <- FetchData(data.combined, vars = GOI,slot = "data") |> 
  rownames_to_column("barcodes") |> 
  pivot_longer(names_to = "gene",values_to = "exp",-barcodes) |> 
  # try to min/max normalize the count varaible per gene in order to rescale the difference in terms of expression
  group_by(gene) %>%
  # threshold of positiveness is based on the distriubtion of the expression of the signal in tihs case
  mutate(exp_min_max = ((exp - min(exp))/(max(exp)-min(exp))),
         exp_cat = case_when(exp > 0~"pos",
                             T~"neg")) %>%
  ungroup() %>%
  mutate(exp_fix = exp + rnorm(nrow(.))/100000)

# get the coordinates
UMAP1_df <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column(var = "barcodes")

# generate the dataset for mapping the data in the umamp
dim(UMAP1_df)
dim(df_exp)
dim(meta)

df_tot <- reduce(list(meta,UMAP1_df,df_exp),left_join, by="barcodes")
df_tot_avg <- df_tot %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined$group <- paste0(data.combined$pathology_class,"-",data.combined$expertAnno.l1,"-",data.combined$orig.ident)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI,group.by = c("group"))

# plot general UMAP -------------------------------------------------------
# build the plot using both info
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = expertAnno.l1),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident.pdf",width = 7,height = 5)

# no lab
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_type2),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_noLab.pdf",width = 7,height = 5)

# split the sample
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = expertAnno.l1),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw() +
  facet_wrap(~pathology_class)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  )
ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_split.pdf",width = 17,height = 15)

# nolab
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_type2),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw() +
  facet_wrap(~pathology_class)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  )
ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_split_noLab.pdf",width = 17,height = 15)

# expression distribution -------------------------------------------------
# crop the 0 expressing cells
df_exp %>%
  ggplot(aes(x=exp))+geom_histogram()+facet_grid(~gene)+theme_bw()+scale_x_log10()+geom_vline(xintercept = 1,col="red",linetype="dotted")

# keep the 0 expressing cells
df_exp %>%
  ggplot(aes(x=exp))+geom_histogram()+facet_wrap(~gene)+theme_bw()+
  # scale_x_log10()+
  geom_vline(xintercept = 2.5,col="red",linetype="dotted")

# library(scales)
# show_col(c("#4662D7FF","#FABA39FF","#7A0403FF"))

# plotting expression -----------------------------------------------------
# by counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~pathology_class) +
  theme_bw() +
  scale_color_gradient(low = "gray",high = "blue") +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_TSPO_count.pdf",width = 16,height = 15)

# by min max normalized counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~pathology_class) + theme_bw() + scale_color_gradient(low = "gray",high = "blue") +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_TSPO_minmax.pdf",width = 16,height = 15)

# plot the category. being 0 or non zero per cell
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.2) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  facet_rep_wrap(gene~pathology_class,repeat.tick.labels = "all",nrow=3)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_bw() + scale_color_manual(values = c("gray","blue")) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_TSPO_propPos.pdf",width = 16,height = 15)

# violin plot fro PTX3 expression use macro categories
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  # this is the processing shown in the violinplot function
  # mutate(exp_fix = exp + rnorm(nrow(.))/100000) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=exp_fix)) + 
  geom_violin(scale = "width")+
  geom_point(position=position_jitter(width = 0.2),alpha=0.1) +
  facet_wrap(~expertAnno.l1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/ManualClean/violin_38_annotationConfident_TSPO.pdf",width = 15,height = 15)

# try to depict the average expression there is roughly one sample per condition
df_avg <- average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  mutate(expertAnno.l1 = str_remove_all(group,pattern = "CX_Ctrl.|CX_Demye.|CX_Mye.|WM_CA.|WM_CI.|WM_Core.|WM_Ctrl.|WM_NAWM.|.s\\d+")) %>%
  # mutate(annotation_confident = str_remove_all(group,pattern = "CX_Ctrl.|CX_Demye.|CX_Mye.|WM_CA.|WM_CI.|WM_Core.|WM_Ctrl.|WM_NAWM.|.s\\d+")) %>%
  mutate(pathology_class = str_extract(group,pattern = "CX_Ctrl|CX_Demye|CX_Mye|WM_CA|WM_CI|WM_Core|WM_Ctrl|WM_NAWM")) %>%
  mutate(sample = str_extract(group,pattern = "s\\d+"))

df_avg |>
  mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core"))) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~expertAnno.l1,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/ManualClean/dotplot_38_annotationConfident_TSPO_expressionAvg.pdf",width = 9,height = 9)

# check is it is a marker in the cell type markes
# library(SeuratWrappers)
# data.combined$group2 <- paste0(data.combined$pathology_class,"-",data.combined$annotation_confident)
# DefaultAssay(data.combined) <- "RNA"

Idents(data.combined) <- "expertAnno.l1"
# subset the object to include only the one of interest

list_test <- lapply(names(table(data.combined$expertAnno.l1)),function(cell_id){
  # check the progress
  print(cell_id)
  test_obj <- subset(data.combined,subset = expertAnno.l1 == cell_id)
  Idents(test_obj) <- "pathology_class"
  test <- RunPrestoAll(test_obj, min.pct = 0.05, logfc.threshold = 0.1)
  
  test |> 
    mutate(cell_id = cell_id)
})

# save a table
df_test_subset <- list_test |> 
  bind_rows()

df_test_subset |> 
  write_tsv("../../out/table/ManualClean/FindAllMarkers_annotation_confident_pathology_class_minpct5_logfcthr01.tsv")

# check if TSPO is in there
df_test_subset |> 
  filter(gene %in% GOI)

# make the stats using a pairwise test on the average expression
# in this case the origin and the sample are exclusives, therefore there is no need to generate another aggregation
df_avg |> 
  finalfit(formula = avg_exp~expertAnno.l1+pathology_class)

# cell_id <- "Astrocytes"
list_pairwise <- lapply(names(table(df_avg$expertAnno.l1)),function(cell_id){
  # check the progress
  print(cell_id)
  # subset the dataset to the single cell tyep
  test <- df_avg |> 
    filter(expertAnno.l1 %in% cell_id)
  
  # run the pairwise test
  test_out <- pairwise.t.test(test$avg_exp, test$pathology_class,p.adjust.method = "fdr")
  
  # wrangle the table
  test2 <- test_out$p.value |> 
    data.frame() |> 
    rownames_to_column("first") |> 
    pivot_longer(names_to = "second",values_to = "padj",-first) |> 
    dplyr::filter(!is.na(padj)) |> 
    mutate(cell_id = cell_id)
  return(test2)
})

list_pairwise |> 
  bind_rows() |>
  arrange(padj)

# # use the crossing option
# crossing(first = data.combined$annotation_confident,second = data.combined$annotation_confident) |> 
#   # remove the equal ones
#   filter(first != second)

# barplot over time of the positive cells for Ptx3 per cell_type
# use the full panel of cells
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  group_by(pathology_class,expertAnno.l1) %>%
  summarise(cells = n(),
            pos = sum(exp_cat=="pos")) %>%
  mutate(prop_pos = pos/cells) %>%
  mutate(log_number_cells = log10(cells)) |> 
  mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core"))) |> 
  ggplot(aes(x=pathology_class,y=prop_pos,fill=log_number_cells))+geom_col()+facet_wrap(~expertAnno.l1)+ theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/ManualClean/barplot_38_annotationConfident_TSPO_propPosCells.pdf",width = 9,height = 9)

# proportion on total sample
# use the full panel of cells
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  group_by(pathology_class,expertAnno.l1) %>%
  summarise(cells = n(),
            pos = sum(exp_cat=="pos")) %>%
  mutate(tot_cells = sum(cells)) %>%
  mutate(prop_pos = pos/tot_cells) %>%
  mutate(log_number_cells = log10(cells)) |> 
  ggplot(aes(x=pathology_class,y=prop_pos,fill=log_number_cells))+geom_col()+facet_wrap(~expertAnno.l1)+ theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/ManualClean/barplot_38_annotationConfident_TSPO_SamplepropPosCells.pdf",width = 9,height = 9)

# explore expression ------------------------------------------------------
data.combined$group_cellType <- paste0(data.combined$origin,"-",data.combined$expertAnno.l1,"-",data.combined$orig.ident)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
# Idents(data.combined) <- "group_cellType"
DefaultAssay(data.combined) <- "RNA"

average_group_cellType <- AverageExpression(data.combined,features = GOI,group.by = c("group_cellType"))

df_avg_cellType <- average_group_cellType$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  mutate(expertAnno.l1 = str_remove_all(group,pattern = "cortex.|wm.|.s\\d+")) %>%
  mutate(origin = str_extract(group,pattern = "cortex|wm")) %>%
  mutate(sample = str_extract(group,pattern = "s\\d+"))

# save the table
df_avg_cellType %>%
  write_tsv("../../out/table/df_avg_cellType_TSPO.tsv")

df_avg_cellType |>
  mutate(log1p = log1p(avg_exp)) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp,col=origin))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/ManualClean/dotplot_38_annotationConfident_TSPO_expressionAvg_cellType.pdf",width = 6,height = 5)

# reorder by median expression and drop the origin variable
df_avg_cellType |>
  mutate(log1p = log1p(avg_exp)) |> 
  group_by(expertAnno.l1) %>%
  mutate(avg_group_exp = median(avg_exp)) %>%
  ungroup() %>%
  mutate(expertAnno.l1 = fct_reorder(expertAnno.l1,avg_group_exp,.desc = T)) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_cowplot()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank()) + scale_y_sqrt()
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/ManualClean/dotplot_38_annotationConfident_TSPO_expressionAvg_cellType_sofi.pdf",width = 6,height = 5)

# reorder the expression by the median expression
df_avg_cellType |>
  mutate(log1p = log1p(avg_exp)) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp,col=origin))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# check is it is a marker in the cell type markes
# library(SeuratWrappers)
Idents(data.combined) <- "expertAnno.l1"
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)

sobj_total_h.markers |>
  write_tsv("../../out/table/ManualClean/FindAllMarkers_annotation_confident_PosOnly_minpct5_logfcthr01.tsv")

sobj_total_h.markers |> 
  filter(gene %in% GOI)

# make the stats using a linear model on the average expression
# in this case the origin and the sample are exclusives, therefore there is no need to generate another aggregation
df_avg_cellType |> 
  finalfit(formula = avg_exp~expertAnno.l1)

df_avg_cellType |>
  data.frame() |> 
  ff_plot(dependent = "avg_exp",explanatory = "expertAnno.l1")

# remove the outliers samples
df_avg_cellType |>
  mutate(log1p = log1p(avg_exp)) |> 
  filter(!(sample %in% c("s4","s26"))) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp,col=origin))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/ManualClean/dotplot_38_annotationConfident_TSPO_expressionAvg_cellType.pdf",width = 6,height = 5)

df_avg_cellType |> 
  filter(!(sample %in% c("s4","s26"))) |> 
  finalfit(formula = avg_exp~expertAnno.l1)

df_avg_cellType |>
  filter(!(sample %in% c("s4","s26"))) |> 
  data.frame() |> 
  ff_plot(dependent = "avg_exp",explanatory = "expertAnno.l1")

# remove the outlier sample in the vascular cells
df_avg_cellType |>
  mutate(log1p = log1p(avg_exp)) |> 
  filter(!(sample %in% c("s4","s26","s22"))) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp,col=origin,label = sample))+
  # geom_text_repel()+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/ManualClean/dotplot_38_annotationConfident_TSPO_expressionAvg_cellType_removeSample.pdf",width = 6,height = 5)

df_avg_cellType |> 
  filter(!(sample %in% c("s4","s26","s22"))) |>
  finalfit(formula = avg_exp~expertAnno.l1)

df_avg_cellType |>
  filter(!(sample %in% c("s4","s26","s22"))) |>
  data.frame() |> 
  ff_plot(dependent = "avg_exp",explanatory = "expertAnno.l1")

# try additive model with disease condition, controls we disease keep cortex and wm together.
# in disease put control vs all. do not separate by location, keep cx and wm together.
LUT_slide <- df_tot |> 
  group_by(orig.ident,pathology_class) |> 
  summarise() |> 
  mutate(disease = case_when(str_detect(pathology_class,pattern = "Ctrl")~"Ctrl",
                             T~"MS"))
df_avg_cellType |> 
  left_join(LUT_slide,by = c("sample"="orig.ident")) |> 
  filter(!(sample %in% c("s4","s26","s22"))) |>
  finalfit(formula = avg_exp~expertAnno.l1+disease)

df_avg_cellType |> 
  left_join(LUT_slide,by = c("sample"="orig.ident")) |> 
  filter(!(sample %in% c("s4","s26","s22"))) |>
  lm(formula = avg_exp~expertAnno.l1+disease) |> 
  summary()

df_avg_cellType |>
  data.frame() |> 
  left_join(LUT_slide,by = c("sample"="orig.ident")) |> 
  filter(!(sample %in% c("s4","s26","s22"))) |>
  ff_plot(dependent = "avg_exp",explanatory = c("expertAnno.l1","disease"))

# plot the data
df_avg_cellType |> 
  left_join(LUT_slide,by = c("sample"="orig.ident")) |> 
  filter(!(sample %in% c("s4","s26","s22"))) |> 
  ggplot(aes(x=expertAnno.l1,y=avg_exp,col=disease))+
  # geom_text_repel()+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/ManualClean/dotplot_38_annotationConfident_TSPO_expressionAvg_cellTypeDisease_removeSample.pdf",width = 6,height = 5)

df_avg_cellType |> 
  left_join(LUT_slide,by = c("sample"="orig.ident")) |> 
  filter(!(sample %in% c("s4","s26","s22"))) |>
  # filter(expertAnno.l1 %in% c("ASTRO","VAS")) |> 
  lm(formula = avg_exp~expertAnno.l1*disease) |> 
  summary()

# try to run  for each single model
lapply(unique(df_avg_cellType$expertAnno.l1), function(x){
  df_avg_cellType |> 
    left_join(LUT_slide,by = c("sample"="orig.ident")) |> 
    filter(!(sample %in% c("s4","s26","s22"))) |>
    filter(expertAnno.l1 %in% x) |> 
    lm(formula = avg_exp~disease) |> 
    tidy() |> 
    mutate(exportAnno.l1 = x)
}) |> 
  bind_rows()

# disease stage -----------------------------------------------------------
data.combined$group_disease <- paste0(data.combined$pathology_class,"-",data.combined$orig.ident)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
# Idents(data.combined) <- "group_cellType"
DefaultAssay(data.combined) <- "RNA"

average_group_disease <- AverageExpression(data.combined,features = GOI,group.by = c("group_disease"))

df_avg_disease <- average_group_disease$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # mutate(annotation_confident = str_remove_all(group,pattern = "CX_Ctrl.|CX_Demye.|CX_Mye.|WM_CA.|WM_CI.|WM_Core.|WM_Ctrl.|WM_NAWM.|.s\\d+")) %>%
  mutate(pathology_class = str_extract(group,pattern = "CX_Ctrl|CX_Demye|CX_Mye|WM_CA|WM_CI|WM_Core|WM_Ctrl|WM_NAWM")) %>%
  mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core"))) |> 
  mutate(sample = str_extract(group,pattern = "s\\d+"))

# save the table of average expression
df_avg_disease |> 
  write_tsv("../../out/table/df_avg_disease.tsv")

df_avg_disease |>
  mutate(log1p = log1p(avg_exp)) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=avg_exp))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/dotplot_38_annotationConfident_TSPO_expressionAvg_disease.pdf",width = 6,height = 5)

# check is it is a marker in the cell type markes
# library(SeuratWrappers)
Idents(data.combined) <- "pathology_class"
sobj_total_h.markers2 <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)

sobj_total_h.markers2 |> 
  filter(gene %in% GOI)

# make the stats using a linear model on the average expression
# in this case the origin and the sample are exclusives, therefore there is no need to generate another aggregation
df_avg_disease |> 
  finalfit(formula = avg_exp~pathology_class)

df_avg_disease |>
  data.frame() |> 
  ff_plot(dependent = "avg_exp",explanatory = "pathology_class")



# -------------------------------------------------------------------------
# martina asked to remove some samples from the overall sample analysis
df_avg_disease |>
  filter(!(sample %in% c("s4","s26"))) |> 
  mutate(log1p = log1p(avg_exp)) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=avg_exp))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/dotplot_38_annotationConfident_TSPO_expressionAvg_disease_sampleRemoved.pdf",width = 6,height = 5)


# check is it is a marker in the cell type markes
# make the stats using a linear model on the average expression
# in this case the origin and the sample are exclusives, therefore there is no need to generate another aggregation
df_avg_disease |> 
  filter(!(sample %in% c("s4","s26"))) |> 
  finalfit(formula = avg_exp~pathology_class)

df_avg_disease |>
  filter(!(sample %in% c("s4","s26"))) |> 
  data.frame() |> 
  mutate(pathology_class2 = factor(pathology_class,levels = c("WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core","CX_Ctrl","CX_Mye","CX_Demye"))) |> 
  ff_plot(dependent = "avg_exp",explanatory = "pathology_class2")


# test pos TSPO vs neg TSPO MG --------------------------------------------
# subset only the MG cells. define the positive and negative cells per TSPO and run a de analysis to compare pull the genes
# subset onlyt he MG cells
sobj_MG <- subset(data.combined,subset = expertAnno.l1 == 'MG')

sobj_MG@meta.data

# for each cell fetch TSPO expression
sobj_MG_TSPO <- FetchData(sobj_MG,"TSPO",slot = "data") |> 
  rownames_to_column("barcodes")

# add the TSPO expression to the original meta
sobj_MG$TSPO <- sobj_MG@meta.data |> 
  rownames_to_column("barcodes") |> 
  left_join(sobj_MG_TSPO) |> 
  pull("TSPO")
  
# build a category form the expression value
sobj_MG$TSPO_cat <- case_when(sobj_MG$TSPO == 0~"neg",
                              T~"pos")

# summarise the number of cells
table(sobj_MG$TSPO_cat)

ggplot(aes(x=TSPO))+geom_histogram()

# run the DE
Idents(sobj_MG) <- "TSPO_cat"

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group 
res_test_TSPO_MG <- RunPresto(object = sobj_MG,ident.1 = "pos",ident.2 = "neg")

# save the table of top markers
res_test_TSPO_MG |> 
  rownames_to_column("gene") |> 
  mutate(cell_id = "MG") |> 
  write_tsv("../../out/table/res_test_TSPO_MG.tsv")

# plot volcano
volcano_tot <- read_tsv("../../out/table/res_test_TSPO_MG.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no"))

ggplot() +
  geom_point(data = volcano_tot%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+theme_bw()+
  ggrepel::geom_text_repel(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene))+theme_bw()+theme(strip.background = element_blank())
ggsave("../../out/image/volcano_test_TSPO_MG.pdf",width = 12,height = 9)

# try a quick and dirty enrichR on the posirive and negatives
# library(tidyverse)
# library(enrichR)
# library(scales)
# library(patchwork)

# run enrichr with the list of genes in the module
# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Cell"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","HDSigDB_Human_2021","Azimuth_Cell_Types_2021","GO_Biological_Process_2023","Descartes_Cell_Types_and_Tissue_2021","CellMarker_Augmented_2021")

# query -------------------------------------------------------------------
# seelct only the clusters with more than 10 genes as degs

# pull the gene names dividing the up regulated from the downregulated
list_genes <- list(list_UP = volcano_tot %>% filter(DE_cat=="up") %>% pull(gene),
                   list_DOWN = volcano_tot %>% filter(DE_cat=="down") %>% pull(gene))

# define the background
# background <- df_modules$feature

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_enrichr <- lapply(list_genes,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_enrichr %>%
  write_tsv("../../out/table/enrichR_test_TSPO_MG.tsv")

# list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")

plot_list_UP <- list_enrichr %>%
  split(f = .$comparison)

# library(scales)
list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_UP)
ggsave("../../out/image/enrichR_test_TSPO_MG.pdf",width = 13,height = 25,limitsize = FALSE)

# test pos TSPO vs neg TSPO VAS -------------------------------------------
# do the same as above but for the vas cells.
# subset only the VAS cells. define the positive and negative cells per TSPO and run a de analysis to compare pull the genes
# subset onlyt he VAS cells
sobj_VAS <- subset(data.combined,subset = expertAnno.l1 == 'VAS')

sobj_VAS@meta.data
dim(sobj_VAS)

# for each cell fetch TSPO expression
sobj_VAS_TSPO <- FetchData(sobj_VAS,"TSPO",slot = "data") |> 
  rownames_to_column("barcodes")

# add the TSPO expression to the original meta
sobj_VAS$TSPO <- sobj_VAS@meta.data |> 
  rownames_to_column("barcodes") |> 
  left_join(sobj_VAS_TSPO) |> 
  pull("TSPO")

# build a category form the expression value
sobj_VAS$TSPO_cat <- case_when(sobj_VAS$TSPO == 0~"neg",
                              T~"pos")

# summarise the number of cells
table(sobj_VAS$TSPO_cat)

sobj_VAS@meta.data |>
  ggplot(aes(x=TSPO))+geom_histogram()

# run the DE
Idents(sobj_VAS) <- "TSPO_cat"

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group 
res_test_TSPO_VAS <- RunPresto(object = sobj_VAS,ident.1 = "pos",ident.2 = "neg")

# save the table of top markers
res_test_TSPO_VAS |> 
  rownames_to_column("gene") |> 
  mutate(cell_id = "VAS") |> 
  write_tsv("../../out/table/res_test_TSPO_VAS.tsv")

# plot volcano
volcano_tot_VAS <- read_tsv("../../out/table/res_test_TSPO_VAS.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no"))

ggplot() +
  geom_point(data = volcano_tot_VAS%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot_VAS%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+theme_bw()+
  ggrepel::geom_text_repel(data = volcano_tot_VAS%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene))+theme_bw()+theme(strip.background = element_blank())
ggsave("../../out/image/volcano_test_TSPO_VAS.pdf",width = 12,height = 9)

# sofia shared a list of genes to highlight
GOI <- c("HLA-A","B2M","IFI27","HLA-C","FTL","HLA-B","BST2","TMEM59","ICAM2","ACTB","ACTG1","FTH1","IGFBP7","PECAM1","VEGFC","ADAMTS1","CD99","CXCL1","CXCL12")
ggplot() +
  geom_point(data = volcano_tot_VAS%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot_VAS%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+theme_bw()+
  ggrepel::geom_text_repel(data = volcano_tot_VAS%>%filter(gene %in% GOI),nudge_x = 0.5,nudge_y = 0.5,min.segment.length = 0,segment.alpha=0.3,aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene,col=DE_cat))+
  theme_bw()+
  theme(strip.background = element_blank())+scale_color_manual(values = c("black","red"))+
  theme(legend.position = "none")
ggsave("../../out/image/volcano_test_TSPO_VAS2.pdf",width = 9,height = 9)

# try a quick and dirty enrichR on the posirive and negatives
# library(tidyverse)
# library(enrichR)
# library(scales)
# library(patchwork)

# run enrichr with the list of genes in the module
# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Cell"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","HDSigDB_Human_2021","Azimuth_Cell_Types_2021","GO_Biological_Process_2023","Descartes_Cell_Types_and_Tissue_2021","CellMarker_Augmented_2021")

# query -------------------------------------------------------------------
# seelct only the clusters with more than 10 genes as degs

# pull the gene names dividing the up regulated from the downregulated
# there are no down genes in this case
list_genes_VAS <- list(list_UP = volcano_tot_VAS %>% filter(DE_cat=="up") %>% pull(gene),
                   list_DOWN = volcano_tot_VAS %>% filter(DE_cat=="down") %>% pull(gene))

# define the background
# background <- df_modules$feature

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_enrichr_VAS <- lapply(list_genes_VAS,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_enrichr_VAS %>%
  write_tsv("../../out/table/enrichR_test_TSPO_VAS.tsv")

# list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")

plot_list_VAS <- list_enrichr_VAS %>%
  split(f = .$comparison)

# library(scales)
list_plot_VAS <- pmap(list(plot_list_VAS,names(plot_list_VAS)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_VAS)
ggsave("../../out/image/enrichR_test_TSPO_VAS.pdf",width = 6,height = 25,limitsize = FALSE)

