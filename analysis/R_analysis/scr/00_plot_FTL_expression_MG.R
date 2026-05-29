# AIM ---------------------------------------------------------------------
# explore the expression of TSPO and FTL in the MG dataset form WM+CX

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
GOI <- c("FTL","TSPO")

table(data.combined@meta.data$cluster_martina)

# check markers -----------------------------------------------------------
# the default thresholds for the markers are too high to have TSPO, therefore reduce them
Idents(data.combined) <- "cluster_martina"
test.markers <- SeuratWrappers::RunPrestoAll(data.combined,
                                             only.pos = F,
                                             min.pct = 0,
                                             logfc.threshold = 0,
                                             return.thresh = 1)

# save the table
test.markers %>%
  write_tsv("../../out/table/00_MG_subcluster_HarmonySample_martinaCluster_markers.tsv")

test.markers %>%
  filter(gene %in% GOI) %>%
  arrange(gene,p_val_adj,cluster)

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
table(data.combined@meta.data$pathology)
table(data.combined@meta.data$sample)

# how many cells are positive per cluster per gene
df_tot %>%
  group_by(cluster_martina,gene,exp_cat) %>%
  summarise(n = n(),.groups = "drop") %>%
  group_by(cluster_martina,gene) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  arrange(gene,cluster_martina,prop) %>%
  filter(exp_cat == "pos") %>%
  print(n = 100)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined$group <- paste0(data.combined$pathology,"-clu_",data.combined$cluster_martina,"-",data.combined$orig.ident,"-",data.combined$origin)
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
  facet_wrap(gene~pathology) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_DPP3_count.pdf",width = 13,height = 12)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,split.by = "pathology",raster = T,order = T,ncol = 3)
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
ggsave("../../out/image/00_density_plot_FTL_TSPO_MG.pdf",width = 8,height = 3)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.05) +
  facet_wrap(gene~pathology) +
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
  facet_wrap(gene~pathology) +
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
  facet_wrap(gene~pathology)+
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
  mutate(pathology_class = str_extract(group,pattern = c("chronic.active|lesion.edge|chronic.inactive|lesion.edge|control.cortex|control.WM|demyelinated.cortex|lesion.core|myelinated.cortex|NAWM"))) |> 
  mutate(donor = str_extract(group,pattern = c("s\\d+"))) |> 
  mutate(cluster_martina = str_extract(group,pattern = c(".clu_.*.")) %>% str_remove_all(pattern = c("s\\d+|clu_")) %>% str_remove_all(pattern = "\\.wm|\\.cortex") %>% str_remove_all(pattern = "\\.")) %>%
  mutate(origin = str_extract(group,pattern = c(".wm|.cortex")) %>% str_remove_all(pattern = c("\\.")))

# check the expression
df_avg %>%
  filter(avg_exp > 10)

# confirm the annotations
df_avg$cluster_martina %>% table()
df_avg$donor %>% table()
df_avg$origin %>% table()
df_avg$pathology_class %>% table()

# are there any outliers?
df_avg %>%
  ggplot(aes(x=1,y = avg_exp)) + geom_boxplot() + facet_wrap(~gene,scales="free")

# is there any NA in the calssification
df_avg$pathology_class %>% is.na() %>% sum()

# plot the average expresison by cell annotation
df_avg |>
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

# split by origin
df_avg |>
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
  facet_wrap(origin~gene,scales = "free")+
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
  facet_wrap(gene~cluster_martina,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/00_dotplot_annotationConfident_DPP3_expressionAvg_treat.pdf",width = 9,height = 9)

# -------------------------------------------------------------------------
# martina was also interested in explorign the relative correlation of expression for the two genes.
# try to check a general correlation across all the averaged samples
df_avg_wide <- df_avg %>%
  pivot_wider(names_from = gene,values_from = avg_exp)

# save the expression pairs values per subcluster
df_avg_wide %>%
  write_tsv("../../out/table/00_df_avg_wide_TSPO_FTL_MG.tsv")

df_avg_wide %>%
  ggplot(aes(x=FTL,y=TSPO)) +
  geom_smooth(method = "lm") +
  geom_point() +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  theme_bw() +
  facet_wrap(~cluster_martina,
             scales="free") +
  theme(strip.background = element_blank())

df_avg_wide %>%
  ggplot(aes(x=FTL,y=TSPO)) +
  geom_smooth(method = "lm") +
  geom_point() +
  # scale_x_continuous(trans = "log1p") +
  # scale_y_continuous(trans = "log1p") +
  theme_bw() +
  facet_wrap(~cluster_martina,
             scales="free") +
  theme(strip.background = element_blank())

df_avg_wide %>%
  ggplot(aes(x=FTL,y=TSPO)) +
  geom_smooth(method = "lm") +
  geom_point() +
  # scale_x_continuous(trans = "log1p") +
  # scale_y_continuous(trans = "log1p") +
  theme_bw() +
  facet_wrap(origin~cluster_martina,
             scales="free",ncol=7) +
  theme(strip.background = element_blank())

df_cor_stats <- df_avg %>%
  split(f = .$cluster_martina) %>%
  lapply(function(x){
    test <- x %>%
      pivot_wider(names_from = gene,values_from = avg_exp)
    
    cor.test(test$FTL,test$TSPO) %>%
      broom::tidy()
  }) %>%
  bind_rows(.id = "cluster_martina")
