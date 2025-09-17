# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(tidyverse)
library(presto)
library(future)
library(ggrepel)
library(enrichR)
library(scales)
library(patchwork)


# read in the data --------------------------------------------------------
# above I did the generation of the object
so <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegration_AllSoupX_test.rds")

# add tha manual annotation from Martina
id_cells <- so@meta.data %>%
  mutate(cell_id = case_when(seurat_clusters %in% c(0,4,14)~"OLIGO",
                           seurat_clusters %in% c(9)~"OPC",
                           seurat_clusters %in% c(3,12)~"ASTRO",
                           seurat_clusters %in% c(5)~"IMMUNE",
                           seurat_clusters %in% c(13)~"LYM",
                           seurat_clusters %in% c(11)~"VAS",
                           seurat_clusters %in% c(1, 2, 10,6)~"EXC NEU",
                           seurat_clusters %in% c(7,8)~"INH NEU")) %>% 
  pull(cell_id)

# add the cell id to the dataset
so$cell_id <- id_cells

# plot the cells with the new manual annotation
DimPlot(so, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'cell_id',raster = T)
ggsave("../../out/image/ManualClean/UMAP_data.combined_WM_CX_harmonySkipIntegAllSoupX_ManualAnnotation.pdf",width = 10,height = 7)

DimPlot(so, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'cell_id',split.by = "origin",raster = T)
ggsave("../../out/image/ManualClean/UMAP_data.combined_WM_CX_harmonySkipIntegAllSoupX_ManualAnnotationSplit.pdf",width = 15,height = 7)

# Identify conserved cell type markers ------------------------------------
# data
DefaultAssay(so) <- "RNA"
so@meta.data

Idents(so)<-paste0(so$pathology_class,".",so$cell_id)

# test --------------------------------------------------------------------
# fro this specific dataset we are also interested into exploring the difference between cluster 0 and cluster 4
Idents(so)

table(so$pathology_class,so$cell_id)

unique(so$cell_id)
unique(so$pathology_class)

list_res_CX_ctrl_ref <- lapply(str_subset(unique(so$pathology_class),pattern = "CX_Ctrl",negate = T),function(pato){
  lapply(unique(so$cell_id),function(cell){
    ident_1 <- paste0(pato,".",cell)
    ident_2 <- paste0("CX_Ctrl.",cell)
    res <- RunPresto(so,
                           ident.1 = ident_1,
                           ident.2 = ident_2,
                           verbose = T,
                           logfc.threshold = 0,only.pos = F)
    return(res)
  })
})
saveRDS(list_res_CX_ctrl_ref,"../../out/object/list_res_CX_ctrl_ref.rds")


list_res_WM_ctrl_ref <- lapply(str_subset(unique(so$pathology_class),pattern = "CX_Ctrl",negate = T),function(pato){
  lapply(unique(so$cell_id),function(cell){
    ident_1 <- paste0(pato,".",cell)
    ident_2 <- paste0("WM_Ctrl.",cell)
    res <- RunPresto(so,
                     ident.1 = ident_1,
                     ident.2 = ident_2,
                     verbose = T,
                     logfc.threshold = 0,only.pos = F)
    return(res)
  })
})
saveRDS(list_res_WM_ctrl_ref,"../../out/object/list_res_WM_ctrl_ref.rds")


# calculate also the average expression per cluster gene
# focus on the genes in the table only
id_genes2 <- rownames(DE_4_vs_0) %>% 
  unique()

avg_exp_cluster2 <- AverageExpression(so,group.by = "seurat_clusters",features = id_genes2)
df_avg_exp_cluster2 <- avg_exp_cluster2$RNA %>%
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(names_to = "cluster",values_to = "avg_exp",-gene) %>% 
  mutate(cluster = str_remove_all(cluster,pattern = "^X")) %>% 
  filter(cluster %in% c("0","4")) %>% 
  pivot_wider(names_from = cluster,values_from = avg_exp,names_prefix = "avg_clust_")

DE_4_vs_0 %>%
  rownames_to_column("gene") %>% 
  left_join(df_avg_exp_cluster2,by=c("gene")) %>% 
  write_tsv("../../out/table/ManualClean/DE_4_vs_0_WM_CX_harmonySkipIntegAllSoupX.tsv")

# plot the volcano
# render all of them as a volcano plot
test_significant <- DE_4_vs_0 %>%
  rownames_to_column("gene") %>% 
  # mutate(comparison = str_extract(annotation,pattern = "MS_vs_CTRL")) %>%
  mutate(threshold = case_when(abs(avg_log2FC) > 1 & p_val_adj<0.05~1,
                               T~0)) %>%
  filter(threshold == 1)

# library(ggrepel)
DE_4_vs_0 %>%
  rownames_to_column("gene") %>% 
  # mutate(comparison = str_extract(annotation,pattern = "MS_vs_CTRL")) %>%
  # filter(symbol %in% setdiff(GOI_SLC,GOI)) %>%
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj))) +
  geom_point(alpha = 0.01) +
  geom_point(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj)),col="red",alpha = 0.5) +
  geom_vline(xintercept = c(-1,1),col="gray",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="gray",linetype="dashed")+
  geom_text_repel(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj),label = gene)) +
  # facet_grid(comparison~cluster) +
  theme_bw()+theme(strip.background = element_blank())
ggsave("../../out/image/ManualClean/volcano_MS_vs_CTRL_manualClean_harmony.pdf",width = 20,height = 5)

# plot a dotplot for all the significant genes in the whole dataset
test_long2 <- DotPlot(so, features = test_significant$gene, dot.scale = 8,cluster.idents = T)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# theme(axis.text.x = element_text(angle = 45,hjust = 1))
test_long2
ggsave("../../out/image/ManualClean/Dotplot_data.combined_WM_CX_harmonySkipIntegAllSoupX.pdf",width = 23,height = 4)

# -------------------------------------------------------------------------
# try to run an enrichment analysis on the top hits
dbs <- listEnrichrDbs()
#filter fo the db of interest
# dbs %>%
#   filter(str_detect(libraryName,pattern = "Atlas"))
# 
# dbs %>%
#   filter(str_detect(libraryName,pattern = "Tabula"))
# 
# dbs %>%
#   filter(str_detect(libraryName,pattern = "Azim"))

dbs_db <- c("KEGG_2021_Human",
            "MSigDB_Hallmark_2020",
            "Reactome_2016",
            "HDSigDB_Human_2021","Tabula_Sapiens","Azimuth_Cell_Types_2021")

# GENE SELECTION ----------------------------------------------------------
list_res_tot_UP <- list("DE_4_vs_0_up" = DE_4_vs_0 %>%
                          rownames_to_column("gene") %>%
                          filter(p_val_adj<0.05,avg_log2FC>1) %>% 
                          pull(gene))

# list_res_tot_DOWN <- list("DE_4_vs_0_down" = DE_4_vs_0 %>%
#                           rownames_to_column("gene") %>%
#                           filter(p_val_adj<0.05,avg_log2FC<(-1)) %>% 
#                           pull(gene))


# query -------------------------------------------------------------------
# seelct only the clusters with more than 10 genes as degs
list_res_tot_UP_filter <- list_res_tot_UP[lengths(list_res_tot_UP)>5]

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_UP <- lapply(list_res_tot_UP_filter,function(x){
  genes <- x
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_UP %>%
  write_tsv("../../out/table/ManualClean/enrichR_test_4_vs_0_UP.tsv")

# explore a specific term
list_UP %>%
  filter(str_detect(Term,pattern = "Olig")) %>% 
  dplyr::select(annotation,Term,Genes)

# list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")

plot_list_UP <- list_UP %>%
  split(f = .$comparison)

# library(scales)
list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
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
ggsave("../../out/image/ManualClean/enrichR_test_4_vs_0_UP.pdf",width = 6,height = 12,limitsize = FALSE)

