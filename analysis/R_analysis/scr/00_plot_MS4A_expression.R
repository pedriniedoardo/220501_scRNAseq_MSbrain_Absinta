# AIM ---------------------------------------------------------------------
# explore the expression of MS4A in the CX+WM dataset, in particular the IMM subset

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
data.combined <- readRDS("../../out/object/129_MG_subcluster_HarmonySample_martinaCluster.rds")
DimPlot(data.combined,label = T,raster = F,group.by = "cluster_martina")

# str_subset(rownames(data.combined),pattern = "HIF")
# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("MS4A3","MS4A6A","MS4A4E","MS4A4A","MS4A6E","MS4A14","MS4A5","MS4A1","MS4A13","MS4A8","MS4A18","MS4A15")

table(data.combined@meta.data$cluster_martina)

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
df_tot_avg <- df_tot %>% group_by(cluster_martina) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)

head(data.combined@meta.data)
table(data.combined@meta.data$pathology_class)
table(data.combined@meta.data$orig.ident)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined$group <- paste0(data.combined$pathology_class,"|",data.combined$cluster_martina,"|",data.combined$orig.ident)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI,group.by = c("group"))

# plot general UMAP -------------------------------------------------------
# build the plot using both info
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=cluster_martina),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cluster_martina),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident.pdf",width = 7,height = 5)

# no lab
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=cluster_martina),size=0.3) +
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
ggsave("../../out/image/00_UMAPggplot_annotationConfident_MS4A_count2.pdf",width = 13,height = 9)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,raster = T,order = T)
# ggsave("../../out/image/06_UMAPSeurat_annotationConfident_DPP3_count.pdf",width = 6,height = 5)

# try nebulosa option
plot_density(data.combined, GOI,reduction = "umap")
ggsave("../../out/image/00_density_plot_MS4A_IMM.pdf",width = 13,height = 9)

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
  ggplot(aes(x=pathology,y=exp_fix)) + 
  geom_violin(scale = "width")+
  geom_point(position=position_jitter(width = 0.2),alpha=0.01) +
  facet_wrap(~cluster_martina) +
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
  mutate(cluster_martina = str_extract(group,pattern = c("cl14|cl17|cl9Absinta|cont|dendritic|homeo.MG|MIMS.foamy|MIMS.iron|mye.MG|per.mac|plasmablast|stress.MG|T.cells|und")))
# there seems to be one outlier sample in the dataset with a very high average experssion
# filter(group != "b_NAWM.clu_1.S20")

df_avg %>%
  filter(avg_exp > 10)

# plot the average expresison by cell annotation
df_avg |>
  # reorder the gene to have a better plot
  group_by(gene) %>%
  mutate(test = mean(avg_exp)) %>%
  ungroup() %>%
  mutate(gene = fct_reorder(gene,test,.desc = T)) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=cluster_martina,y=avg_exp))+
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
ggsave("../../out/image/00_boxplot_cluster_AvgExp_MS4A_IMM.pdf",width = 12,height = 8)

# do the same as above but split by condition
df_avg |>
  # reorder the gene to have a better plot
  group_by(gene) %>%
  mutate(test = mean(avg_exp)) %>%
  ungroup() %>%
  mutate(gene = fct_reorder(gene,test,.desc = T)) %>%
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
ggsave("../../out/image/00_boxplot_treat_AvgExp_MS4A_IMM.pdf",width = 8,height = 10)

# rank the cluster_martina category by the average expression per sample (regardelss of the gene)
rank_cell_id <- df_avg %>%
  group_by(cluster_martina) %>%
  summarise(avg = mean(avg_exp)) %>%
  arrange(desc(avg)) %>%
  pull(cluster_martina)

# generate the average estimates per cell annotation per gene
summary_exp_cell_id <- df_avg %>%
  group_by(cluster_martina, gene) %>%
  summarise(avg = mean(avg_exp)) %>%
  mutate(cluster_martina = factor(cluster_martina,levels = rank_cell_id))

# plot splitting by treat full
df_avg |>
  mutate(cluster_martina = factor(cluster_martina,levels = rank_cell_id)) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=pathology_class,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_grid(gene~cluster_martina,scales = "free") +
  geom_hline(data = summary_exp_cell_id,aes(yintercept = avg),col="red",linetype = "dashed")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/00_boxplot_TreatCluster_AvgExp_MS4A_IMM.pdf",width = 15,height = 12)

# dotplot -----------------------------------------------------------------
# make a list out of the genes
GOI_list <- list(GOI = GOI)

test_dotplot <- DotPlot(data.combined,
                        features = unique(unlist(GOI_list)),
                        dot.scale = 8,
                        cluster.idents = T,
                        group.by = "cluster_martina") +
  RotatedAxis() +
  labs(title = "cluster_martina") +
  theme(strip.text = element_text(angle = 90))

df_test <- lapply(GOI_list,function(x){
  test_dotplot$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

test_long01 <- df_test %>%
  # force the order
  # mutate(id = factor(id,levels = c(9,5,6,3,2,7,4,8,0,1,10))) %>% 
  # mutate(cell_type = factor(cell_type,levels = c("Bc","Tc","RBc","FIBRO","MAC","DC","ENDO"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 0))+
  scale_color_gradient(low = "lightgrey",high = "blue",limits = c(-1,2),oob = scales::squish)
ggsave(plot = test_long01,"../../out/image/00_dotplot_cluster_MS4A_IMM.pdf",width = 6,height = 5)
