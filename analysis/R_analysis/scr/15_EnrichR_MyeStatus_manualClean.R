# libraries ---------------------------------------------------------------
library(tidyverse)
library(enrichR)
library(scales)
library(patchwork)

# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Azim"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","HDSigDB_Human_2021","GTEx_Aging_Signatures_2021","ClinVar_2019","GO_Biological_Process_2023")

# GENE SELECTION ----------------------------------------------------------
# list of genes to consider for the enrichment analysis
# test_out <- readRDS("out/IMMUNE_CCA_martina_trimm_KO_C1QB.rds")

# read in the metadata for the genes
# meta <- read_csv("data/scTrem/GSE130626_gene_info.csv")

# # save the ranked object, also change the genes into genenames
# file_id <- dir("out/table/") %>%
#   str_subset(pattern = "FindAllMarkers_data.combined_fix_") %>%
#   str_subset(pattern = "top100",negate = T)
#
# name <- case_when(file_id=="FindAllMarkers_data.combined_fix_DoubletSinglet_seuratClusters.tsv"~"basic",
#                   file_id=="FindAllMarkers_data.combined_fix_regressCC_DoubletSinglet_seuratClusters.tsv"~"regressCC",
#                   file_id=="FindAllMarkers_data.combined_fix_regressCCRIBO_DoubletSinglet_seuratClusters.tsv"~"regressCCRIBO" )

folder <- "../../out/table/ManualClean/"
file <- dir(folder) %>%
  str_subset(pattern = "responsePseudobulkShr_MyeStatus")

df_res <- lapply(file, function(x){
  test_plot <- read_tsv(paste0(folder,x))
}) %>%
  setNames(file) %>%
  bind_rows(.id = "file") %>%
  mutate(treat = str_extract(file,pattern = "DeMye_vs_Norm|ReMye_vs_Norm|DeMye_vs_ReMye"),
         cell_id = str_extract(file,pattern = "clust_0|clust_2|clust_5|clust_6|clust_13")) %>%
  mutate(comparison = paste0(treat,"|",cell_id))

list_df <- split(df_res,f = df_res$comparison)

list_res_tot_UP <- lapply(list_df, function(x){
  x %>%
    filter(padj<0.05,log2FoldChange>1) %>%
    pull(gene)
})

list_res_tot_DOWN <- lapply(list_df, function(x){
  x %>%
    filter(padj<0.05,log2FoldChange<1) %>%
    pull(gene)
})

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
  write_tsv("../../out/table/ManualClean/enrichR_Pseudobulk_MyeStatus_UP.tsv")

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
ggsave("../../out/image/ManualClean/enrichR_Pseudobulk_MyeStatus_UP.pdf",width = 20,height = 40,limitsize = FALSE)

# query -------------------------------------------------------------------
# seelct only the clusters with more than 10 genes as degs
list_res_tot_DOWN_filter <- list_res_tot_DOWN[lengths(list_res_tot_DOWN)>5]

list_DOWN <- lapply(list_res_tot_DOWN_filter,function(x){
  genes <- x
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_DOWN %>%
  write_tsv("../../out/table/ManualClean/enrichR_Pseudobulk_MyeStatus_DOWN.tsv")

# list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")

plot_list_DOWN <- list_DOWN %>%
  split(f = .$comparison)

# library(scales)
list_plot_DOWN <- pmap(list(plot_list_DOWN,names(plot_list_DOWN)), function(x,y){
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

wrap_plots(list_plot_DOWN)
ggsave("../../out/image/ManualClean/enrichR_Pseudobulk_MyeStatus_DOWN.pdf",width = 20,height = 40,limitsize = FALSE)

