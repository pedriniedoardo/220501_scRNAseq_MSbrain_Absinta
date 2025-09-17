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
library(enrichR)
library(scales)
library(patchwork)

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
# mono.de <- FindMarkers(sobj, ident.1 = "EXC NEU|CX_Demye", ident.2 = "EXC NEU|CX_Ctrl", verbose = T)
df_comparison_sc <- bind_rows(
  data.frame(ident1 = paste0("EXC NEU|",c("CX_Demye","CX_Mye","CX_Demye")),
             ident2 = paste0("EXC NEU|",c("CX_Ctrl","CX_Ctrl","CX_Mye"))) ,
  
  data.frame(ident1 = paste0("INH NEU|",c("CX_Demye","CX_Mye","CX_Demye")),
             ident2 = paste0("INH NEU|",c("CX_Ctrl","CX_Ctrl","CX_Mye"))) 
)

sc.de <- pmap(df_comparison_sc,function(ident1,ident2){
  print(paste(ident1,"vs",ident2))
  
  mono.de <- SeuratWrappers::RunPresto(sobj, ident.1 = ident1, ident.2 = ident2, verbose = T,logfc.threshold = 0)
  df_final <- mono.de %>%
    rownames_to_column("gene") %>%
    mutate(ident1 = ident1,
           ident2 = ident2,
           compar = paste(ident1,"vs",ident2))
  
  return(df_final)
}) %>%
  bind_rows()

table(sc.de$compar)

# head(mono.de)
# The p-values obtained from this analysis should be interpreted with caution, because these tests treat each cell as an independent replicate and ignore inherent correlations between cells originating from the same sample. Such analyses have been shown to find a large number of false positive associations, as has been demonstrated by Squair et al., 2021, Zimmerman et al., 2021, Junttila et al., 2022, and others. Below, we show how pseudobulking can be used to account for such within-sample correlation.

# save the table
sc.de %>%
  write_tsv("../../out/table/140_sc.de.tsv")

# sample plot
sc.de %>%
  mutate(cell_id = str_extract(ident1,"EXC NEU|INH NEU")) %>%
  tibble() %>%
  mutate(test = str_extract_all(compar,"CX_Demye|CX_Ctrl|CX_Mye")) %>%
  mutate(test2 = map(test,function(x){paste0(x,collapse = " vs ")}) %>% unlist()) %>%
  select(-test) %>%
  mutate(test_sig = case_when(abs(avg_log2FC)>1 & p_val_adj < 0.05 ~ "sig",
                              T~"not sig")) %>%
  ggplot(aes(x=avg_log2FC,y=-log(p_val_adj))) +
  geom_point(shape = 1,alpha=0.4,aes(col = test_sig)) +
  facet_grid(test2~cell_id) +
  geom_vline(xintercept = c(-1,1),col="red",linetype = "dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype = "dashed")+
  theme_bw() +
  scale_color_manual(values = c("black","red")) +
  theme(strip.background = element_blank())

# alternative plot reducing the threshold of significance
sc.de %>%
  mutate(cell_id = str_extract(ident1,"EXC NEU|INH NEU")) %>%
  tibble() %>%
  mutate(test = str_extract_all(compar,"CX_Demye|CX_Ctrl|CX_Mye")) %>%
  mutate(test2 = map(test,function(x){paste0(x,collapse = " vs ")}) %>% unlist()) %>%
  select(-test) %>%
  mutate(test_sig = case_when(abs(avg_log2FC)>0.5 & p_val_adj < 0.05 ~ "sig",
                              T~"not sig")) %>%
  ggplot(aes(x=avg_log2FC,y=-log(p_val_adj))) +
  geom_point(shape = 1,alpha=0.4,aes(col = test_sig)) +
  facet_grid(test2~cell_id) +
  geom_vline(xintercept = c(-0.5,0.5),col="red",linetype = "dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype = "dashed")+
  theme_bw() +
  scale_color_manual(values = c("black","red")) +
  theme(strip.background = element_blank())

# run enrichR on the genes that have abs FC > 0.5 -------------------------

# run enrichr with the list of genes in the module
# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
# filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "GO"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_Pathways_2024","GO_Biological_Process_2023")

# query -------------------------------------------------------------------

list_genes_UP <- read_tsv("../../out/table/140_sc.de.tsv") %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5 & !is.na(gene)) %>%
  split(f = .$compar) %>%
  map(function(x){
    x %>%
      pull(gene)
  })

list_genes_DOWN <- read_tsv("../../out/table/140_sc.de.tsv") %>%
  filter(p_val_adj < 0.05 & avg_log2FC < -(0.5) & !is.na(gene)) %>%
  split(f = .$compar) %>%
  map(function(x){
    x %>%
      pull(gene)
  })

# # pull the gene names dividing the up regulated from the downregulated
# list_genes <- list(list_UP = results$res_GMPvsHUVEC_shr %>% filter(log2FoldChange>0) %>% pull(symbol),
#                    list_DOWN = results$res_GMPvsHUVEC_shr %>% filter(log2FoldChange<0) %>% pull(symbol))

# define the background
# background <- df_modules$feature

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_enrichr_UP <- lapply(list_genes_UP,function(x){
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

list_enrichr_UP %>%
  write_tsv("../../out/table/140_enrichR_scDE_NEU_UP.tsv")

list_enrichr_DOWN <- lapply(list_genes_DOWN,function(x){
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

list_enrichr_DOWN %>%
  write_tsv("../../out/table/140_enrichR_scDE_NEU_DOWN.tsv")

# list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")
plot_list_UP <- list_enrichr_UP %>%
  split(f = .$comparison)

plot_list_DOWN <- list_enrichr_DOWN %>%
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

wrap_plots(list_plot_UP,nrow = 1)
ggsave("../../out/image/140_enrichR_scDE_NEU_UP.pdf",width = 36,height = 12,limitsize = FALSE)

list_plot_DOWN <- pmap(list(plot_list_DOWN,names(plot_list_DOWN)), function(x,y){
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

wrap_plots(list_plot_DOWN,nrow = 1)
ggsave("../../out/image/140_enrichR_scDE_NEU_DOWN.pdf",width = 36,height = 12,limitsize = FALSE)

# explore the DE tables ---------------------------------------------------
# Martina asked to provide the numbers of genes that are differentially expressed across dataset
read_tsv("../../out/table/140_sc.de.tsv") %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5 & !is.na(gene)) %>%
  group_by(compar) %>%
  summarise(n = n())
  
read_tsv("../../out/table/140_sc.de.tsv") %>%
  filter(p_val_adj < 0.05 & avg_log2FC < (-0.5) & !is.na(gene)) %>%
  group_by(compar) %>%
  summarise(n = n())
