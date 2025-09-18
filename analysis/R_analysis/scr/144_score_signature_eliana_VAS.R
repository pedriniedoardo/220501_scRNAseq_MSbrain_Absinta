# AIM ---------------------------------------------------------------------
# this is the updated script for the generation of the modules scores for eliana's custom signatures
# I am running the analysis on the VAS subset only 

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)
library(msigdbr)
library(UpSetR)

# read the data -----------------------------------------------------------
# read in the dataset
data.combined <- readRDS("../../out/object/130_VAS_subcluster_HarmonySample.rds")

# confirm the object identity
DimPlot(data.combined,raster=T,group.by = "expertAnno.l1")
DimPlot(data.combined,raster=F,group.by = "RNA_snn_res.0.4",label = T)

# load the siganture file
list_sig <- readRDS("../../out/object/143_custom_signatures_eliana_endo.rds")

# wrangling ---------------------------------------------------------------
# add the new classification to the metadata
meta <- data.combined@meta.data %>%
  rownames_to_column("barcodes")

# meta_full <- left_join(meta,LUT,by=c("official_id"))
meta_full <- meta

# add to the original dataset
# data.combined$pathology_class <- meta_full$pathology_class

# score the siganture in the UMAP -----------------------------------------
# run the enrichment for the signature. do it on the UMAP using the score siganatures
DefaultAssay(data.combined) <- "RNA"
# x <- "ECFC_MSvsCTRL_POS"

# run the snippet over the whole signatures
list_data2 <- lapply(names(list_sig),function(x){
  # extract the dataframe of genes in the signature
  signature.genes.df <- list_sig[[x]]
  
  # pull the genes
  signature.genes <- signature.genes.df %>%
    pull(Genes) %>%
    unique()
  
  # score the module
  data.combined <- AddModuleScore(data.combined,
                                  features = list(signature.genes),
                                  name="signature_score")
  
  # confirm the addition of the score for the module
  # data.combined@meta.data
  df_meta <- data.combined@meta.data %>%
    rownames_to_column("barcode") %>% 
    mutate(signature = x) 
  # mutate(pathology_class = factor(pathology_class,levels = c("control cortex","myelinated cortex","demyelinated cortex")))
  
  # save the table with the scores
  df_meta %>% 
    write_tsv(paste0("../../out/table/modules_custom_VAS/144_Module_score_",x,".tsv"))
  
  # save the UMAP coordinates
  df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
    data.frame() %>%
    rownames_to_column("barcode")
  
  # data2 <- left_join(df_UMAP,df_meta,"barcode")
  # data2_avg <- data2 %>% group_by(RNA_snn_res.0.4) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
  data2 <- left_join(df_UMAP,df_meta,"barcode") %>%
    mutate(test = paste0(origin,"_",disease))
  data2_avg <- data2 %>% group_by(RNA_snn_res.0.4) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = RNA_snn_res.0.4)) +
    facet_grid(~test) + theme_bw() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/modules_custom_VAS/144_UMAP_score_",x,".pdf"),width = 13,height = 4)
  ggsave(paste0("../../out/image/modules_custom_VAS/144_UMAP_score_",x,".png"),width = 13,height = 4)
  
  # plot the score as distribution but as ridges
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    ggplot(aes(x=signature_score1,y=test,fill=test))+
    ggridges::geom_density_ridges(alpha=0.5)+
    facet_wrap(~RNA_snn_res.0.4,scales = "free")+
    # facet_grid(~paste0(origin,"_",disease)) + theme_bw() + 
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))
  # scale_fill_manual(values = c("green","yellow","red"))
  ggsave(paste0("../../out/image/modules_custom_VAS/144_dist_score_ridges_",x,".pdf"),width = 13,height = 12)
  
  # plot the distribution of the score per cluster
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot(aes(x=test,y = signature_score1,fill=test)) + 
    geom_violin() +
    geom_boxplot(width=0.1,outlier.shape = NA) +
    facet_wrap(~RNA_snn_res.0.4) + theme_bw() + 
    # facet_wrap(pathology_class~orig.ident) + theme_bw() + 
    # scale_fill_manual(values = c("green","yellow","red"))+
    theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))
  ggsave(paste0("../../out/image/modules_custom_VAS/144_violin_score_",x,".pdf"),width = 13,height = 12)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CA","WM_CI","WM_Core"))) %>% 
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    # geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = RNA_snn_res.0.4)) +
    facet_wrap(~test,nrow = 2) + theme_void() + 
    theme(strip.background = element_blank()) +
    scale_color_gradientn("sig score",
                          colours = c("gray","orange","red"),
                          oob = scales::squish,limits = c(0.05,0.25))
  # scale_color_gradientn("sig score",
  #                       colours = c("gray","orange","red"))
  # scale_color_gradientn(colours = c("gray","yellow","red"))
  # scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/modules_custom_VAS/144_UMAP_score_",x,"_tailored1.pdf"),width = 9,height = 8)
  ggsave(paste0("../../out/image/modules_custom_VAS/144_UMAP_score_",x,"_tailored1.png"),width = 9,height = 8,bg="white")
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CA","WM_CI","WM_Core"))) %>% 
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    # geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = RNA_snn_res.0.4)) +
    facet_wrap(~test,nrow = 3) + theme_void() + 
    theme(strip.background = element_blank()) +
    scale_color_gradientn("sig score",
                          colours = viridis::turbo(10),limits = c(0.05,0.2),oob = scales::squish)
  ggsave(paste0("../../out/image/modules_custom_VAS/144_UMAP_score_",x,"_tailored2.pdf"),width = 9,height = 8)
  ggsave(paste0("../../out/image/modules_custom_VAS/144_UMAP_score_",x,"_tailored2.png"),width = 9,height = 8,bg="white")
  
  return(data2)
})

# custom plotting ---------------------------------------------------------
# plot the UMAP for all the signatures in one plot
# x<-list_data2[[1]]
list_plot <- lapply(list_data2,function(x){
  
  # pick the bottom 30% and to max value and use it as threshold for filtering the color scale
  thr_low <- quantile(x$signature_score1, prob=0.30)
  thr_high <- max(x$signature_score1)
  
  thr_low
  thr_high
  
  x %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = RNA_snn_res.0.4)) +
    facet_grid(signature~test) + theme_cowplot() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn("sig score",
                          colours = viridis::turbo(10),limits = c(thr_low,thr_high),oob = scales::squish)
})

wrap_plots(list_plot,ncol=1)
ggsave("../../out/image/144_UMAP_custom_eliana_VAS.pdf",height = 16,width = 17)
ggsave("../../out/image/144_UMAP_custom_eliana_VAS.png",width = 17,height = 16,bg="white")

# load the module score for all the custom signatures
file_sig <- dir("../../out/table/modules_custom_VAS/")

df_modules <- lapply(file_sig,function(x){
  df_modules <- read_tsv(paste0("../../out/table/modules_custom_VAS/",x)) %>%
    mutate(test = paste0(origin,"_",disease))
  return(df_modules)
}) %>%
  bind_rows()


# define the threshold per cell type for the signature
df_90 <- df_modules %>% 
  # filter(signature %in% "ECFC_MSvsCTRL_POS",
  #        test %in% "cortex_CTRL") %>%
  filter(test %in% "cortex_CTRL") %>%
  mutate(RNA_snn_res.0.4 = factor(RNA_snn_res.0.4)) %>%
  group_by(signature, RNA_snn_res.0.4) %>% 
  summarise(thr = quantile(signature_score1, prob=0.90))

# show the distributino of the score per cell type and condition
df_modules %>%
  # focus only on MG and include only baseline Fe and myelin 
  filter(RNA_snn_res.0.4 %in% c("0","5","9")) %>%
  # mutate(BraakStage=as.factor(BraakStage)) %>% 
  ggplot(aes(x=signature_score1,y=test))+
  ggridges::geom_density_ridges(alpha=0.5)+
  facet_wrap(signature~RNA_snn_res.0.4,scales = "free_x",ncol=3)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) +
  # add the threshold for the specific cell type
  geom_vline(data = df_90 %>% filter(RNA_snn_res.0.4 %in% c("0","5","9")), aes(xintercept=thr),linetype = "dashed",col="red")
# scale_fill_manual(values = c("green","yellow","red"))
ggsave(paste0("../../out/image/144_dist_score_ridges_eliana_VASsubclusters_tailored.pdf"),width = 8,height =8)
