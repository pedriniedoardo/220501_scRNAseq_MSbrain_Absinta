# AIM ---------------------------------------------------------------------
# explore the expression of TSPO in the WM + CX dataset.
# replicate the same processing used for CHIT1

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(cowplot)
library(patchwork)
library(Nebulosa)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1",split.by = "origin")
DimPlot(data.combined,label = T,raster = T,group.by = "origin")
DimPlot(data.combined,label = T,raster = T,group.by = "pathology_class")

# str_subset(rownames(data.combined),pattern = "HIF")
# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("FTL","TSPO")

table(data.combined@meta.data$expertAnno.l1)

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

df_tot <- purrr::reduce(list(meta,UMAP1_df,df_exp),left_join, by="barcodes")
df_tot_avg <- df_tot %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)


head(data.combined@meta.data)

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
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident.pdf",width = 7,height = 5)

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
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_noLab.pdf",width = 7,height = 5)

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
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_count.pdf",width = 13,height = 12)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,split.by = "pathology_class",raster = T,order = T,ncol = 3)
# ggsave("../../out/image/06_UMAPSeurat_annotationConfident_DPP3_count.pdf",width = 25,height = 3)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(~gene) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_count2.pdf",width = 6,height = 5)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,raster = T,order = T)
# ggsave("../../out/image/06_UMAPSeurat_annotationConfident_DPP3_count.pdf",width = 6,height = 5)

# try nebulosa option
plot_density(data.combined, GOI,reduction = "umap")

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.05) +
  facet_wrap(gene~pathology_class) +
  theme_cowplot() +
  # scale_color_gradient(low = "gray",high = "blue") +
  scale_color_viridis_c(option = "turbo") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_count_alt.pdf",width = 13,height = 12)

# by min max normalized counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~pathology_class) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_minmax.pdf",width = 13,height = 12)

# plot the category. being 0 or non zero per cell
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.05) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  # facet_rep_wrap(gene~treat,repeat.tick.labels = "all",nrow=3)+
  facet_wrap(gene~pathology_class)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_cowplot() +
  scale_color_manual(values = c("gray","blue")) +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_proppos.pdf",width = 13,height = 12)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.05) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  # facet_rep_wrap(gene~treat,repeat.tick.labels = "all",nrow=3)+
  facet_wrap(~gene)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_cowplot() +
  scale_color_manual(values = c("gray","blue")) +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_proppos2.pdf",width = 6,height = 5)

# violin plot for GOI expression use macro categories
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  # this is the processing shown in the violinplot function
  # mutate(exp_fix = exp + rnorm(nrow(.))/100000) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=exp_fix)) + 
  geom_violin(scale = "width")+
  geom_point(position=position_jitter(width = 0.2),alpha=0.01) +
  facet_wrap(~expertAnno.l1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("../../out/image/00_violin_annotationConfident_DPP3.pdf",width = 15,height = 10)

# try to depict the average expression there is roughly one sample per condition
df_avg <- average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(pathology_class = str_extract(group,pattern = c("CX_Ctrl|CX_Demye|CX_Mye|WM_CA|WM_CI|WM_Core|WM_Ctrl|WM_NAWM"))) |> 
  mutate(donor = str_extract(group,pattern = c("s\\d+"))) |> 
  mutate(expertAnno.l1 = str_extract(group,pattern = c("AST|EPENDYMA|EXC|NEU|IMM|INH|NEU|LYM|OLIGO|OPC|VAS")))

# plot the average expresison by cell annotation
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free")+
  scale_y_continuous(trans = "log1p")

# do the same as above but split by condition
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free") +
  scale_y_continuous(trans = "log1p")

# rank the expertAnno.l1 category by the average expression per sample (regardelss of the gene)
rank_cell_id <- df_avg %>%
  group_by(expertAnno.l1) %>%
  summarise(avg = mean(avg_exp)) %>%
  arrange(desc(avg)) %>%
  pull(expertAnno.l1)

# generate the average estimates per cell annotation per gene
summary_exp_cell_id <- df_avg %>%
  group_by(expertAnno.l1, gene) %>%
  summarise(avg = mean(avg_exp)) %>%
  mutate(expertAnno.l1 = factor(expertAnno.l1,levels = rank_cell_id))

# plot splitting by treat full
df_avg |>
  mutate(expertAnno.l1 = factor(expertAnno.l1,levels = rank_cell_id)) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_grid(gene~expertAnno.l1,scales = "free") +
  geom_hline(data = summary_exp_cell_id,aes(yintercept = avg),col="red",linetype = "dashed")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/06_dotplot_annotationConfident_DPP3_expressionAvg_treatFull.pdf",width = 9,height = 9)

# plot splitting by treat
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(gene~expertAnno.l1,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/00_dotplot_annotationConfident_DPP3_expressionAvg_treat.pdf",width = 9,height = 9)

# -------------------------------------------------------------------------
# martina was also interested in explorign the relative correlation of expression for the two genes.
# try to check a general correlation across all the averaged samples
df_avg_wide <- df_avg %>%
  pivot_wider(names_from = gene,values_from = avg_exp)

df_avg_wide %>%
  ggplot(aes(x=FTL,y=TSPO)) +
  geom_smooth(method = "lm") +
  geom_point() +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  theme_bw() +
  facet_wrap(~expertAnno.l1,
             scales="free") +
  theme(strip.background = element_blank())

df_cor_stats <- df_avg %>%
  split(f = .$expertAnno.l1) %>%
  lapply(function(x){
    test <- x %>%
      pivot_wider(names_from = gene,values_from = avg_exp)
    
    cor.test(test$FTL,test$TSPO) %>%
      broom::tidy()
  }) %>%
  bind_rows(.id = "expertAnno.l1")

# -------------------------------------------------------------------------
# explore the dataset only for the ENDO subset
data.combined_endo <- readRDS("../../out/object/130_VAS_subcluster_HarmonySample.rds")

DimPlot(data.combined_endo,group.by = "RNA_snn_res.0.4",label = T)
plot_density(data.combined_endo, GOI_endo,reduction = "umap")

# track also the high confidence endo cells
meta_endo <- data.combined_endo@meta.data %>%
  mutate(endo_confidence = case_when(RNA_snn_res.0.4 %in% c(0,5,9)~"endo",
                                     T~"other")) %>%
  select(endo_confidence)

# add the new meta to the main metadata
data.combined_endo_update <- AddMetaData(data.combined_endo,meta_endo)

DimPlot(data.combined_endo_update,group.by = "endo_confidence",label = T)

# average the expression by disease and region
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined_endo_update$group <- paste0(data.combined_endo_update$origin,"-",
                                          data.combined_endo_update$disease,"-",
                                          data.combined_endo_update$expertAnno.l1,"-",
                                          data.combined_endo_update$orig.ident)
data.combined_endo_update$group2 <- paste0(data.combined_endo_update$origin,"-",
                                           data.combined_endo_update$disease,"-",
                                           data.combined_endo_update$endo_confidence,"-",
                                           data.combined_endo_update$orig.ident)

# for dge
data.combined_endo_update$location_disease <- paste0(data.combined_endo_update$origin,"-",
                                                     data.combined_endo_update$disease)
data.combined_endo_update$location_disease_endo_confidence <- paste0(data.combined_endo_update$origin,"-",
                                                                     data.combined_endo_update$disease,"-",
                                                                     data.combined_endo_update$endo_confidence)
data.combined_endo_update$disease_endo_confidence <- paste0(data.combined_endo_update$disease,"-",
                                                            data.combined_endo_update$endo_confidence)


# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
# Idents(data.combined_endo_update) <- "group"
# DefaultAssay(data.combined_endo_update) <- "RNA"

GOI_endo <- c("TSPO")

average_GOI_endo <- AverageExpression(data.combined_endo_update,features = GOI_endo,group.by = c("group"))
average_GOI_endo2 <- AverageExpression(data.combined_endo_update,features = GOI_endo,group.by = c("group2"))

df_avg_endo <- average_GOI_endo$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI_endo) %>%
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(location_class = str_extract(group,pattern = c("cortex|wm"))) %>% 
  mutate(disease_class = str_extract(group,pattern = c("CTRL|MS"))) %>%
  mutate(donor = str_extract(group,pattern = c("s\\d+"))) %>%
  mutate(disease2 = paste(location_class,disease_class,sep = "|"))

df_avg_endo2 <- average_GOI_endo2$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI_endo) %>%
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(location_class = str_extract(group,pattern = c("cortex|wm"))) %>% 
  mutate(disease_class = str_extract(group,pattern = c("CTRL|MS"))) %>%
  mutate(donor = str_extract(group,pattern = c("s\\d+"))) %>%
  mutate(endo_class = str_extract(group,pattern = c("endo|other"))) %>%
  mutate(disease2 = paste(location_class,disease_class,sep = "|"))

# split by conidtion and location
p01 <- df_avg_endo |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=disease2,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free")+
  scale_y_continuous(trans = "log1p") +
  labs(title = "Endo subset")

p012 <- df_avg_endo2 |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=disease2,y=avg_exp,color = endo_class))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,jitter.height = 0,dodge.width = 0.8),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free")+
  scale_y_continuous(trans = "log1p") +
  labs(title = "Endo subset")


# do the same as above but split by condition
p02 <- df_avg_endo |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=disease_class,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free")+
  scale_y_continuous(trans = "log1p") +
  labs(title = "Endo subset")

p022 <- df_avg_endo2 |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=disease_class,y=avg_exp,color = endo_class))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,jitter.height = 0,dodge.width = 0.8),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free")+
  scale_y_continuous(trans = "log1p") +
  labs(title = "Endo subset")

p02 | p01
p022 | p012

# -------------------------------------------------------------------------
# add the meta from endo subset to the main annotation
# the idea is to check if the expression of TSPO comes mainly from the high confidence endo subset of the other vascular cells
data.combined_test <- AddMetaData(data.combined,meta_endo)

# add the meta to the main object
data.combined_test$location_disease_cellid_endo_confidence <- paste0(data.combined_test$origin,"-",
                                                              data.combined_test$disease,"-",
                                                              data.combined_test$expertAnno.l1,"-",
                                                              data.combined_test$endo_confidence,"-",
                                                              data.combined_test$orig.ident)
# used for DGE at the single cell level
data.combined_test$cellid_endo_confidence <- paste0(data.combined_test$expertAnno.l1,"-",
                                                    data.combined_test$endo_confidence)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
average_GOI_test <- AverageExpression(data.combined_test,features = GOI_endo,group.by = c("location_disease_cellid_endo_confidence"))

df_avg_test <- average_GOI_test$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI_endo) %>%
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(location_class = str_extract(group,pattern = c("cortex|wm"))) %>% 
  mutate(disease_class = str_extract(group,pattern = c("CTRL|MS"))) %>%
  mutate(donor = str_extract(group,pattern = c("s\\d+"))) %>%
  mutate(endo_class = str_extract(group,pattern = c("endo|other"))) %>%
  mutate(expertAnno.l1 = str_extract(group,pattern = c("AST|EPENDYMA|EXC|NEU|IMM|INH|NEU|LYM|OLIGO|OPC|VAS"))) %>%
  mutate(disease2 = paste(location_class,disease_class,sep = "|"))

# save the table of average expression
df_avg_test %>%
  write_tsv("../../out/table/00_avgExp_data.combined_endo_confidence.tsv")

df_avg_test |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp,color = endo_class))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,jitter.height = 0,dodge.width = 0.8),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free")+
  scale_y_continuous(trans = "log1p")

# try markers at the sc level
DefaultAssay(data.combined_test) <- "RNA"
Idents(data.combined_test) <- "cellid_endo_confidence"
table(Idents(data.combined_test))

sobj_total_markers_test <- RunPrestoAll(data.combined_test,
                                        only.pos = F,
                                        min.pct = 0.01,
                                        logfc.threshold = 0)

# check the GOIs
sobj_total_markers_test %>%
  filter(gene %in% GOI_endo)

sobj_total_markers_test %>%
  write_tsv("../../out/table/00_markers_data.combined_endo_confidence.tsv")

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# make the comparison for the subset of VAS cells
# the idea is to identify the main source of the expressio of TSPO from the VAS subcluster
DefaultAssay(data.combined_endo_update) <- "RNA"
Idents(data.combined_endo_update) <- "RNA_snn_res.0.4"
table(Idents(data.combined_endo_update))

(DimPlot(data.combined_endo,group.by = "RNA_snn_res.0.4",label = T) + ggtitle("Harmony Sample"))
(DimPlot(data.combined_endo,split.by = "disease",group.by = "RNA_snn_res.0.4",label = T) + ggtitle("Harmony Sample"))
# (DimPlot(data.combined_endo,split.by = "location_disease") + ggtitle("Harmony Sample"))

# small marker panel to justify the identity of each RNA_snn_res.0.4 cluster (vascular subtypes + likely contaminants picked up from the broader WM+CX dataset)
# 0: capillary/BBB endothelial | 5: venous endothelial | 9: arterial endothelial
# 4: ambiguous/non-specific (ribosomal+metallothionein-high, no coherent canonical EC/pericyte/SMC signal - likely a high-depth or stress/IFN substate)
# 2: pericyte, capillary-type (PDGFRB/ABCC9/CSPG4-high) | 6: pericyte, precapillary/ensheathing-type (RGS5/NDUFA4L2/HIGD1B/MUSTN1-high)
# 7: smooth muscle cell | 1: perivascular fibroblast | 3: neuron | 8: oligodendrocyte | 10: astrocyte | 11: immune/microglia | 12: OPC
marker_panel_VAS <- list(
  "Endo" = c("CLDN5","PECAM1","FLT1"),
  "Art EC" = c("BMX","SEMA3G","EFNB2"),
  "Ven EC" = c("VWF","NR2F2"),
  "Peri" = c("RGS5","KCNJ8","ABCC9","HIGD1B","NDUFA4L2","PDGFRB"),
  "SMC" = c("ACTA2","TAGLN","MYH11"),
  "Fibro" = c("COL1A2","DCN","FBLN1"),
  "OPC" = c("PDGFRA","CSPG4","VCAN","PTPRZ1"),
  "Oligo" = c("PLP1","MOG","MBP"),
  "Astro" = c("AQP4","GFAP","SLC1A2"),
  "Neu" = c("SNAP25","SYT1"),
  "Imm" = c("PTPRC","CSF1R","P2RY12")
)

test <- DotPlot(data.combined_endo_update, features = marker_panel_VAS, group.by = "RNA_snn_res.0.4")
# ggsave("../../out/image/00_dotplot_VAS_subcluster_res0.4_markers.pdf",width = 12,height = 5)

df_test <- lapply(marker_panel_VAS,function(x){
  test$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

test_long01 <- df_test %>%
  # force the order
  mutate(id = factor(id,levels = c(0,9,5,4,2,6,7,1,12,8,10,3,11))) %>% 
  mutate(cell_type = factor(cell_type,levels = c("Endo","Art EC","Ven EC","Peri","SMC","Fibro","OPC","Oligo","Astro","Neu","Imm"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue",limits = c(-1,2),oob = scales::squish)
ggsave(plot = test_long01,filename = "../../out/image/00_dotplot_res0.4_VAS_subcluster.pdf",width = 10,height = 5)

# priduce the markers
sobj_total_h.res0.4 <- RunPrestoAll(data.combined_endo_update,
                                    only.pos = F,
                                    min.pct = 0.01,
                                    logfc.threshold = 0)

# save the table of all markers
sobj_total_h.res0.4 %>%
  write_tsv("../../out/table/00_markers_res0.4_VAS_subcluster.tsv")

# check gois
sobj_total_h.res0.4 %>%
  filter(gene %in% GOI_endo)

# -------------------------------------------------------------------------
# test by disease
DefaultAssay(data.combined_endo_update) <- "RNA"
Idents(data.combined_endo_update) <- "location_disease"
table(Idents(data.combined_endo_update))

# find markers for every cluster compared to all remaining cells, picj the reference and run aver all the other levels
pct2.id_vec <- c("wm-CTRL","cortex-CTRL")

list_res2 <- lapply(pct2.id_vec, function(pct2.id){
  # track the progress
  print(pct2.id)
  
  # build the comparisons
  pct1.id_vec <- table(Idents(data.combined_endo_update)) %>% names() %>% str_subset(negate = T,pattern = pct2.id)
  
  # make the comparisons
  list_res <- lapply(pct1.id_vec,function(pct1.id){
    # track the progress
    print(pct1.id)
    
    # run the DGE
    sobj_total_h.disease <- RunPresto(data.combined_endo_update,
                                      ident.1 = pct1.id,
                                      ident.2 = pct2.id,
                                      only.pos = F,
                                      min.pct = 0.01,
                                      logfc.threshold = 0) %>%
      rownames_to_column("gene") %>%
      mutate(pct.1_id = pct1.id,
             pct.2_id = pct2.id,
             test = paste0(pct.1_id,"_vs_",pct2.id))
  })
  
  # make a data.frame
  list_res %>%
    bind_rows()
  
})

# make a single data.frame
df_res <- list_res2 %>%
  bind_rows()

df_res$test %>% table()

df_res %>%
  filter(gene %in% c("TSPO"))

# save the table of all markers
df_res %>%
  write_tsv("../../out/table/00_DEG_disease_VAS_subcluster_disease.tsv")

# -------------------------------------------------------------------------
# test by disease endo class and location
DefaultAssay(data.combined_endo_update) <- "RNA"
Idents(data.combined_endo_update) <- "location_disease_endo_confidence"
table(Idents(data.combined_endo_update))

# find markers for every cluster compared to all remaining cells, picj the reference and run aver all the other levels
pct2.id_vec2 <- c("wm-CTRL-endo","cortex-CTRL-endo")

list_res3 <- lapply(pct2.id_vec2, function(pct2.id){
  # track the progress
  print(pct2.id)
  
  # build the comparisons
  pct1.id_vec <- table(Idents(data.combined_endo_update)) %>% names() %>% str_subset(negate = T,pattern = pct2.id)
  
  # make the comparisons
  list_res <- lapply(pct1.id_vec,function(pct1.id){
    # track the progress
    print(pct1.id)
    
    # run the DGE
    sobj_total_h.disease <- RunPresto(data.combined_endo_update,
                                      ident.1 = pct1.id,
                                      ident.2 = pct2.id,
                                      only.pos = F,
                                      min.pct = 0.01,
                                      logfc.threshold = 0) %>%
      rownames_to_column("gene") %>%
      mutate(pct.1_id = pct1.id,
             pct.2_id = pct2.id,
             test = paste0(pct.1_id,"_vs_",pct2.id))
  })
  
  # make a data.frame
  list_res %>%
    bind_rows()
  
})

# make a single data.frame
df_res2 <- list_res3 %>%
  bind_rows()

df_res2$test %>% table()

df_res2 %>%
  filter(gene %in% c("TSPO"))

# save the table of all markers
df_res2 %>%
  write_tsv("../../out/table/00_DEG_disease_VAS_subcluster_disease2.tsv")


# -------------------------------------------------------------------------
# test by disease and endo class
DefaultAssay(data.combined_endo_update) <- "RNA"
Idents(data.combined_endo_update) <- "disease_endo_confidence"
table(Idents(data.combined_endo_update))

# find markers for every cluster compared to all remaining cells, picj the reference and run aver all the other levels
pct2.id_vec3 <- c("CTRL-endo")

list_res4 <- lapply(pct2.id_vec3, function(pct2.id){
  # track the progress
  print(pct2.id)
  
  # build the comparisons
  pct1.id_vec <- table(Idents(data.combined_endo_update)) %>% names() %>% str_subset(negate = T,pattern = pct2.id)
  
  # make the comparisons
  list_res <- lapply(pct1.id_vec,function(pct1.id){
    # track the progress
    print(pct1.id)
    
    # run the DGE
    sobj_total_h.disease <- RunPresto(data.combined_endo_update,
                                      ident.1 = pct1.id,
                                      ident.2 = pct2.id,
                                      only.pos = F,
                                      min.pct = 0.01,
                                      logfc.threshold = 0) %>%
      rownames_to_column("gene") %>%
      mutate(pct.1_id = pct1.id,
             pct.2_id = pct2.id,
             test = paste0(pct.1_id,"_vs_",pct2.id))
  })
  
  # make a data.frame
  list_res %>%
    bind_rows()
  
})

# make a single data.frame
df_res3 <- list_res4 %>%
  bind_rows()

df_res3$test %>% table()

df_res2 %>%
  filter(gene %in% c("TSPO"))

# save the table of all markers
df_res2 %>%
  write_tsv("../../out/table/00_DEG_disease_VAS_subcluster_disease2.tsv")
