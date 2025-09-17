# AIM ---------------------------------------------------------------------
# the aim of this test is to compare seurat's implementaiton of pseudobulk analysis, to the default process of DGE
# the reference of the test is presented here: https://satijalab.org/seurat/articles/de_vignette.html
# this second script run the comparison on the DGE

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(tidyverse)
library(harmony)
library(ggExtra)
library(ComplexUpset)
library(cowplot)
library(UpSetR)
library(DESeq2)

# read in the dataset -----------------------------------------------------
sobj <- readRDS("../../out/object/129_MG_subcluster_HarmonySample_martinaCluster.rds")

# check the object version
class(sobj@assays$RNA)

# test DGE at single cell -------------------------------------------------
# add in one covariate the cell anntation and the stimulation status
# in this case I need only the disease status

# set the ident
Idents(sobj) <- "disease"

# run the DGE over the same cell type for stim vs ctrl. The units are the single cells.
# log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
mono.de <- FindMarkers(sobj, ident.1 = "MS", ident.2 = "CTRL", verbose = FALSE)

head(mono.de)
# The p-values obtained from this analysis should be interpreted with caution, because these tests treat each cell as an independent replicate and ignore inherent correlations between cells originating from the same sample. Such analyses have been shown to find a large number of false positive associations, as has been demonstrated by Squair et al., 2021, Zimmerman et al., 2021, Junttila et al., 2022, and others. Below, we show how pseudobulking can be used to account for such within-sample correlation.

# save the table
saveRDS(mono.de,"../../out/object/02_mono.de.rds")

# test DGE at pseudobulk --------------------------------------------------
# pseudobulk the counts based on donor-condition-celltype
pseudo_sobj <- AggregateExpression(sobj, assays = "RNA", return.seurat = T, group.by = c("dataset", "disease"),slot = "counts")

# confirm the count slot is correctly populated, notice that the slot contains integers
pseudo_sobj@assays$RNA@counts
pseudo_sobj@assays$RNA@data

# add the covariate for the stimulation per cell type
# pseudo_sobj$sample.disease <- paste(pseudo_sobj$dataset, pseudo_sobj$disease, sep = "_")
pseudo_sobj$sample.disease <- 
  pseudo_sobj@meta.data %>%
  rownames_to_column("sample.disease") %>%
  pull("sample.disease")

pseudo_sobj$disease <- pseudo_sobj$sample.disease %>%
  str_extract(pattern = "MS|CTRL")

# explore the dimensionality of the new dataset
pseudo_sobj@meta.data %>%
  group_by(disease) %>%
  summarise(n = n())

# run the DGE over the same cell type for MS vs CTRL This time the unist are the pseudobulks per donor/celltype/stimulus
Idents(pseudo_sobj) <- "disease"
bulk.mono.de <- FindMarkers(object = pseudo_sobj, 
                            ident.1 = "MS", 
                            ident.2 = "CTRL",
                            test.use = "DESeq2")
head(bulk.mono.de)

# save the table
saveRDS(bulk.mono.de,"../../out/object/02_bulk.mono.de.rds")

# -------------------------------------------------------------------------
# is there any gene from the Senmayo
table <- read_tsv("../../data/signatures/senescence/senmayo_signature_short_fixed_review.txt")

mono.de %>%
  rownames_to_column("gene") %>%
  filter(gene %in% table$human_gene)

bulk.mono.de %>%
  rownames_to_column("gene") %>%
  filter(gene %in% table$human_gene)

# compare FC --------------------------------------------------------------
# merge all the stat per gene
df_full <- mono.de %>%
  rownames_to_column("gene") %>%
  left_join(bulk.mono.de %>%
              rownames_to_column("gene"),
            by = c("gene"),
            suffix = c(".sc",".pBulk"))

# compare the estimated fc for both analysis
df_full %>%
  ggplot(aes(x=avg_log2FC.sc,y=avg_log2FC.pBulk))+geom_point(shape=1)+theme_bw()+geom_abline(slope = 1,intercept = 0,col="red",linetype = "dashed")
ggsave("../../out/image/02_compare_sc_pbulk_FC.pdf",width = 4,height = 4)

# see the script

# compare p values --------------------------------------------------------
# compare the estimated adjusted p-value for both analysis
p_pval <- df_full %>%
  ggplot(aes(x=p_val_adj.sc,y=p_val_adj.pBulk))+geom_point(shape=1) +
  theme_bw() +
  geom_abline(slope = 1,intercept = 0,col="red",linetype = "dashed") +
  coord_fixed()

ggMarginal(p_pval, type="histogram")

# try with the raw pvalues
df_full %>%
  ggplot(aes(x=p_val.sc,y=p_val.pBulk))+
  # geom_point(shape=1) +
  coord_fixed() +
  # scale_y_continuous(trans = "log1p") +
  # scale_x_continuous(trans = "log1p") +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE,n = 500) +
  # theme_cowplot()+
  # scale_fill_viridis_c(option = "turbo",trans = "log1p")+
  scale_fill_viridis_c(
    trans = "log1p", 
    limits = c(NA, 10),  # Apply the max cutoff
    oob = scales::squish,                        # Squish values exceeding the cutoff
    name = "Log Density",option = "turbo"
  ) +
  geom_abline(slope = 1,intercept = 0,col="red",linetype = "dashed") +
  theme_cowplot() +
  theme(
    axis.line = element_blank()      # Keep axis lines
  )

# try to apply the filters for significance to compare the number of genes
# try the upset plot version 
# make a list of the significnat genes. make sure to devide coherent up and coherent down in both dataset
list_sig_up <- lapply(list(sc = mono.de,pBulk = bulk.mono.de), function(x){
  x %>%
    rownames_to_column("gene") %>%
    filter(p_val_adj < 0.05) %>%
    filter(avg_log2FC>1) %>%
    pull(gene)
})

list_sig_down <- lapply(list(sc = mono.de,pBulk = bulk.mono.de), function(x){
  x %>%
    rownames_to_column("gene") %>%
    filter(p_val_adj < 0.05) %>%
    filter(avg_log2FC< (-1)) %>%
    pull(gene)
})

df_sig_up <- lapply(list(sc = mono.de,pBulk = bulk.mono.de), function(x){
  x %>%
    rownames_to_column("gene") %>%
    filter(p_val_adj < 0.05) %>%
    filter(avg_log2FC>1)
}) %>%
  bind_rows(.id = "test")

# plot option 1
# UpSetR::upset(fromList(list_sig_up), order.by = "freq") 

# plot option 2
(ComplexUpset::upset(fromList(list_sig_up),colnames(fromList(list_sig_up)),wrap=T) + ggtitle("genes up")) +
  (ComplexUpset::upset(fromList(list_sig_down),colnames(fromList(list_sig_down)),wrap=T) + ggtitle("genes down"))

# who is the single gene down and up
list_test <- list(down = list_sig_down,
                  up = list_sig_up)

list_int <- lapply(list_test, function(test){
  
  # pull the elements in the lists
  df1 <- lapply(test,function(x){
    data.frame(gene = x)
  }) %>% 
    bind_rows(.id = "path")
  
  # head(df1)
  
  # pull the inique features
  df2 <- data.frame(gene=unique(unlist(test)))
  
  # head(df2)
  
  # define the intersections for each feature across all the elemenet in the list
  df_int <- lapply(df2$gene,function(x){
    # pull the name of the intersections
    intersection <- df1 %>% 
      dplyr::filter(gene==x) %>% 
      arrange(path) %>% 
      pull("path") %>% 
      paste0(collapse = "|")
    
    # build the dataframe
    data.frame(gene = x,int = intersection)
  }) %>% 
    bind_rows()
  
  # head(df_int,n=20)
  
  return(df_int)
})

# confirm the summaries from upset
lapply(list_int, function(df_int){
  df_int %>% 
    group_by(int) %>% 
    summarise(n=n()) %>% 
    arrange(desc(n))
})


# pull the genes from a specific intersection
lapply(list_int, function(df_int){
  df_int %>%
    filter(int %in% c("pBulk"))
})

# -------------------------------------------------------------------------
# explore some genes that are labelled as DEGs in both versions of the analysis
# top 5
top5_common <- df_full %>%
  arrange(p_val.pBulk) %>%
  filter(abs(avg_log2FC.pBulk)>1) %>%
  filter(abs(avg_log2FC.sc)>1) %>%
  dplyr::slice(1:5) %>%
  pull(gene)

# plot top5 by condition
VlnPlot(sobj, features = top5_common, group.by = "disease",ncol = 5)
