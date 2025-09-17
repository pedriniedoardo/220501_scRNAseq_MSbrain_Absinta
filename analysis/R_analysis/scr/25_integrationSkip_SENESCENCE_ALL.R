# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)

# read the data -----------------------------------------------------------
# read in the dataset
data.combined <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegration_SCtypeAnnotation.rds")

# load the new classification from martina
# LUT <- read_csv("../../data/LUT_sample_MyeClass.csv")

# load the siganture file
list_sig <- readRDS("../../data/signatures/senescence_pathways.rds")

# wrangling ---------------------------------------------------------------
# add the new classification to the metadata
# meta <- data.combined@meta.data %>%
#   rownames_to_column("barcodes")
# 
# meta_full <- left_join(meta,LUT,by=c("official_id"))
# 
# # add to the original dataset
# data.combined$pathology_class <- meta_full$pathology_class

# score the siganture in the UMAP -----------------------------------------
# run the enrichment for the signature. do it on the UMAP using the score siganatures
DefaultAssay(data.combined) <- "RNA"
# x <- "senmayo"

# run the snippet over the whole signatures
lapply(names(list_sig),function(x){
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
    mutate(signature = x) %>% 
    mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","WM_Ctrl","WM_NAWM","CX_Mye","WM_CI","WM_CA","CX_Demye","WM_Core")))
  
  # save the table with the scores
  df_meta %>% 
    write_tsv(paste0("../../out/table/ManualClean/modules_SENESCENCE/Module_score_CX_WM_",x,".tsv"))
  
  # save the UMAP coordinates
  df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
    data.frame() %>%
    rownames_to_column("barcode")
  
  # data2 <- left_join(df_UMAP,df_meta,"barcode")
  # data2_avg <- data2 %>% group_by(seurat_clusters) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
  data2 <- left_join(df_UMAP,df_meta,"barcode")
  data2_avg <- data2 %>% group_by(seurat_clusters) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
  
  # data2 %>%
  #   arrange(signature_score1) %>%
  #   # mutate(gene = "Ptx3") %>%
  #   ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
  #   geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  #   facet_grid(~pathology) + theme_bw() + 
  #   theme(strip.background = element_blank()) +
  #   # scale_color_gradientn(colours = c("blue","gray","red"))
  #   scale_color_gradientn(colours = viridis::turbo(10))
  # # scale_color_gradientn(colours = c("blue","gray","red"), 
  # #                       values = rescale(c(-0.1,0,0.58)),
  # #                       guide = "colorbar", limits=c(-0.1,0.58))
  # ggsave(paste0("out/image/modules/UMAP_score_",x,".pdf"),width = 12,height = 3)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
    facet_grid(~pathology_class) + theme_bw() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/UMAP_score_CX_WM_",x,".pdf"),width = 12,height = 3)
  ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/UMAP_score_CX_WM_",x,".png"),width = 9,height = 3)
  
  # data2 %>%
  #   arrange(signature_score1) %>%
  #   # mutate(gene = "Ptx3") %>%
  #   ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
  #   geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  #   facet_wrap(pathology~patient) + theme_bw() + 
  #   theme(strip.background = element_blank()) +
  #   # scale_color_gradientn(colours = c("blue","gray","red"))
  #   scale_color_gradientn(colours = viridis::turbo(10))
  # ggsave(paste0("out/image/modules/UMAP_score_split_",x,".pdf"),width = 12,height = 9)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
    facet_wrap(pathology_class~orig.ident) + theme_bw() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/UMAP_score_split_CX_WM_",x,".pdf"),width = 29,height = 28)
  ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/UMAP_score_split_CX_WM_",x,".png"),width = 29,height = 28)
  
  # plot the score as distribution
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    ggplot(aes(x=signature_score1,fill=pathology_class))+geom_density(alpha=0.5)+
    facet_wrap(~seurat_clusters,scales = "free")+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))
    # scale_fill_manual(values = c("green","yellow","red"))
  ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/dist_score_CX_WM_",x,".pdf"),width = 12,height =9)
  
  # plot the score as distribution but as ridges
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    ggplot(aes(x=signature_score1,y=pathology_class,fill=pathology_class))+
    ggridges::geom_density_ridges(alpha=0.5)+
    facet_wrap(~seurat_clusters,scales = "free")+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))
    # scale_fill_manual(values = c("green","yellow","red"))
  ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/dist_score_ridges_CX_WM_",x,".pdf"),width = 12,height =9)
  
  # plot the distribution of the score per cluster
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    mutate(name_sample = paste0(pathology_class,"_",orig.ident)) %>% 
    mutate(rank = rank(pathology_class,ties.method = "min")) %>% 
    mutate(name_sample = fct_reorder(name_sample, rank,.desc = F)) %>% 
    
    ggplot(aes(x=name_sample,y = signature_score1,fill=pathology_class)) + 
    geom_violin() +
    geom_boxplot(width=0.1,outlier.shape = NA) +
    facet_wrap(~seurat_clusters) + theme_bw() + 
    # scale_fill_manual(values = c("green","yellow","red"))+
    theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))
  ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/violin_score_split_CX_WM_",x,".pdf"),width = 25,height = 9)
  
  # plot the distribution of the score per cluster
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot(aes(x=pathology_class,y = signature_score1,fill=pathology_class)) + 
    geom_violin() +
    geom_boxplot(width=0.1,outlier.shape = NA) +
    facet_wrap(~seurat_clusters) + theme_bw() + 
    #mscale_fill_manual(values = c("green","yellow","red"))+
    theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))
  ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/violin_score_CX_WM_",x,".pdf"),width = 9,height = 9)
})
