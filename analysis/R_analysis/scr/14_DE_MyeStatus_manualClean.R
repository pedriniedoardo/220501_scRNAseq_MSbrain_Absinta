# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(limma)
library(ggrepel)

# read in the data --------------------------------------------------------
# read in the annotated object
data.combined <- readRDS("../../out/object/ManualClean/data.combined_fix_filter_norm_doublet_integrated_SoupX_manualClean_harmony.rds")

# load the new classification from martina
LUT <- read_csv("../../data/LUT_sample_MyeClass.csv")

# wrangling ---------------------------------------------------------------
# add the new classification to the metadata
meta <- data.combined@meta.data %>%
  rownames_to_column("barcodes")

meta_full <- left_join(meta,LUT,by=c("official_id"))

# add to the original dataset
data.combined$pathology_class <- meta_full$pathology_class

# compare the clusters per treatment --------------------------------------
# make suer the correct defalult dataset il loaded should be RNA
DefaultAssay(data.combined)

# check that in the RNA slot the data object is indeed loaded with normalized values
data.combined@assays$RNA@data[1:10,1:10]

# define the grouping variables for the comparison of CSF vs CTRL
head(data.combined@meta.data)

# define the new grouping
table(data.combined@meta.data$seurat_clusters,data.combined@meta.data$pathology_class)

data.combined$condition.annotation <- paste(data.combined$pathology_class,data.combined$seurat_clusters, sep = "_")
head(data.combined@meta.data)

# update the idents of the object
Idents(data.combined) <- "condition.annotation"

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
# run the differential expression over the WT
# clusters_id <- as.character(sort(unique(data.combined$annotation_confident))) %>%
#   str_subset(pattern = "P09",negate = T)
clusters_id <- as.character(sort(unique(data.combined$seurat_clusters)))

# demyelinated cortex vs control cortex
# run presto to make the comparison faster
list_DeMye_vs_Norm <- lapply(clusters_id,function(x){
  id_1 <- paste0("demyelinated cortex_",x)
  id_2 <- paste0("control cortex_",x)
  response <- RunPresto(data.combined, ident.1 = id_1, ident.2 = id_2, verbose = T,logfc.threshold = 0)
  response %>%
    rownames_to_column("gene") %>%
    mutate(id_1 = id_1,
           id_2 = id_2) %>%
    mutate(cluster = x)
})

list_DeMye_vs_Norm %>%
  setNames(paste0(clusters_id,"_DeMye_vs_Norm")) %>%
  bind_rows(.id = "annotation") %>%
  write_tsv("../../out/table/ManualClean/response_DeMye_vs_Norm_manualClean_harmony.tsv")

# remyelinated cortex vs control cortex
# run presto to make the comparison faster
list_ReMye_vs_Norm <- lapply(clusters_id,function(x){
  id_1 <- paste0("myelinated cortex_",x)
  id_2 <- paste0("control cortex_",x)
  response <- RunPresto(data.combined, ident.1 = id_1, ident.2 = id_2, verbose = T,logfc.threshold = 0)
  response %>%
    rownames_to_column("gene") %>%
    mutate(id_1 = id_1,
           id_2 = id_2) %>%
    mutate(cluster = x)
})

list_ReMye_vs_Norm %>%
  setNames(paste0(clusters_id,"_ReMye_vs_Norm")) %>%
  bind_rows(.id = "annotation") %>%
  write_tsv("../../out/table/ManualClean/response_ReMye_vs_Norm_manualClean_harmony.tsv")

# demyelinated cortex vs control cortex
# run presto to make the comparison faster
list_DeMye_vs_ReMye <- lapply(clusters_id,function(x){
  id_1 <- paste0("demyelinated cortex_",x)
  id_2 <- paste0("myelinated cortex_",x)
  response <- RunPresto(data.combined, ident.1 = id_1, ident.2 = id_2, verbose = T,logfc.threshold = 0)
  response %>%
    rownames_to_column("gene") %>%
    mutate(id_1 = id_1,
           id_2 = id_2) %>%
    mutate(cluster = x)
})

list_DeMye_vs_ReMye %>%
  setNames(paste0(clusters_id,"_DeMye_vs_ReMye")) %>%
  bind_rows(.id = "annotation") %>%
  write_tsv("../../out/table/ManualClean/response_DeMye_vs_ReMye_manualClean_harmony.tsv")

# -------------------------------------------------------------------------
# test <- read_tsv("out/table/response_treat_vs_P09_regressCC_DoubletSinglet.tsv")
# # check the expression of some of the genes on interes
# GOI <- c("IRF7","DDX58","CD34","PROCR")
# test %>%
#   filter(gene %in% GOI)
# plot the volcanos per cluster -------------------------------------------
folder <- "../../out/table/ManualClean/"
file <- dir(folder) %>%
  str_subset(pattern = "response_") %>% 
  str_subset(pattern = "MS_vs_CTRL",negate = T)

df_res <- lapply(file, function(x){
  test_plot <- read_tsv(paste0(folder,x))
}) %>%
  bind_rows()

# show the distribution of the pvalues
df_res %>%
  mutate(comparison = str_extract(annotation,pattern = "DeMye_vs_Norm|ReMye_vs_Norm|DeMye_vs_ReMye")) %>%
  ggplot(aes(x = p_val)) + geom_histogram()+theme_bw()+facet_grid(comparison~cluster)+theme(strip.background = element_blank())
ggsave("../../out/image/ManualClean/dist_p_value_MyeStatus_manualClean_harmony.pdf",width = 20,height = 12)

# render all of them as a volcano plot
test_significant <- df_res %>%
  mutate(comparison = str_extract(annotation,pattern = "DeMye_vs_Norm|ReMye_vs_Norm|DeMye_vs_ReMye")) %>%
  mutate(threshold = case_when(abs(avg_log2FC) > 1 & p_val_adj<0.05~1,
                               T~0)) %>%
  filter(threshold == 1)

# library(ggrepel)
df_res %>%
  mutate(comparison = str_extract(annotation,pattern = "DeMye_vs_Norm|ReMye_vs_Norm|DeMye_vs_ReMye")) %>%
  # filter(symbol %in% setdiff(GOI_SLC,GOI)) %>%
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj))) +
  geom_point(alpha = 0.01) +
  geom_point(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj)),col="red",alpha = 0.5) +
  geom_vline(xintercept = c(-1,1),col="gray",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="gray",linetype="dashed")+
  geom_text_repel(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj),label = gene)) +
  facet_grid(comparison~cluster) +
  theme_bw()+theme(strip.background = element_blank())
ggsave("../../out/image/ManualClean/volcano_MyeStatus_manualClean_harmony.pdf",width = 20,height = 15)
