# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(limma)
library(DESeq2)

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

# confirm the groupign for the metadata
table(data.combined@meta.data$official_id,data.combined@meta.data$seurat_clusters,data.combined@meta.data$disease)
table(data.combined@meta.data$official_id,data.combined@meta.data$seurat_clusters)

# wragling ----------------------------------------------------------------
# create the column for the aggragation
data.combined$pseudobulk <- paste0(data.combined@meta.data$official_id,".",
                                   "clust_",data.combined@meta.data$seurat_clusters,".",
                                   data.combined@meta.data$disease)

# save the aggregtated table of counts
matrix_counts <- AggregateExpression(object = data.combined,
                                     group.by = c("pseudobulk"),
                                     assays = "RNA",
                                     slot = "counts",
                                     return.seurat = FALSE) %>%
  .$RNA

dim(matrix_counts)
matrix_counts[1:10,1:10]

# save the matrix as object
matrix_counts %>%
  saveRDS("../../out/object/ManualClean/pseudobulk_MSStatus_SoupX_manualClean_harmony.rds")

# save the aggregated mete
meta_summary <- data.combined@meta.data %>%
  group_by(official_id,disease,sex,age,pathology_class,orig.ident.cca,pseudobulk) %>%
  summarise()
meta_summary %>%
  write_tsv("../../out/table/ManualClean/pseudobulk_Meta_MSStatus_SoupX_manualClean_harmony.tsv")

# run the sample deseq analysis -------------------------------------------
# run the analysis only on some clusters 0,2,5,6,13
table(data.combined@meta.data$seurat_clusters,data.combined@meta.data$official_id)

# Let's run DE analysis
# 1. Get counts matrix
list_count <- lapply(c("clust_0","clust_2","clust_5","clust_6","clust_13"),function(x){
  matrix_counts %>%
    data.frame() %>% 
    rownames_to_column("gene") %>% 
    select(gene,contains(x)) %>%
    column_to_rownames("gene")
}) %>%
  setNames(c("clust_0","clust_2","clust_5","clust_6","clust_13"))

# counts_astro <- df_counts %>%
#   select(gene,contains("Astrocytes")) %>%
#   column_to_rownames("gene")

# 2. generate sample level metadata
# read in the summarised meta
meta_summary <- read_tsv("../../out/table/ManualClean/pseudobulk_Meta_MSStatus_SoupX_manualClean_harmony.tsv")

list_colData <- lapply(list_count,function(data_table){
  data.frame(samples = colnames(data_table)) %>%
    separate(samples,into = c("id_sample","cell_type","treat1"),sep = "\\.",remove = F) %>%
    # mutate(sample_id_fix = str_sub(string = id_sample,start = 2,end = -1)) %>%
    left_join(meta_summary,by = c("samples"="pseudobulk")) %>%
    column_to_rownames(var = 'samples')
})

# colData <- data.frame(samples = colnames(counts_astro)) %>%
#   separate(samples,into = c("id_sample","cell_type","treat1","treat2"),sep = "\\.",remove = F) %>%
#   mutate(sample_id_fix = str_sub(string = id_sample,start = 2,end = -1)) %>%
#   left_join(meta_summary,by = c("sample_id_fix"="orig.ident.cca")) %>%
#   column_to_rownames(var = 'samples')

# build the design
list_design <- lapply(list_colData,function(x){
  treat <- factor(x$disease,levels = c("CTRL","MS"))
  gender <- factor(x$sex)
  # age <- x$age

  # design <- model.matrix(~ gender+age+treat)
  design <- model.matrix(~ gender+treat)
  # colnames(design) <- c("intercept","genderM","age","treatMS")
  colnames(design) <- c("intercept","genderM","treatMS")
  design
})

# save the disegn
# saveRDS(design,"out/object/design.rds")

# perform DESeq2 --------
# Create DESeq2 object
list_dds <- pmap(list(list_count,list_colData,list_design),function(x,y,z){
  DESeqDataSetFromMatrix(countData = x,
                         colData = y,
                         design = z)
})

# filter
# x <- list_dds$clust_0
# name_dds <- "clust_0"
# design_dds <- list_design$clust_0

pmap(list(list_dds,names(list_dds),list_design),function(x,name_dds,design_dds){
  dds <- x
  keep <- rowSums(counts(dds)) >=10
  dds <- dds[keep,]

  # run DESeq2
  dds2 <- DESeq(dds)

  # Check the coefficients for the comparison
  resultsNames(dds2)

  # save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
  contrast <- makeContrasts(MSvsCTRL = treatMS,
                            levels = design_dds)

  # Generate results object
  res_MS <- results(dds2, contrast=contrast[,"MSvsCTRL"],alpha = 0.05)

  summary(res_MS)

  res_MS %>% data.frame() %>% rownames_to_column("gene") %>% arrange(pvalue) %>%
    write_tsv(paste0("../../out/table/ManualClean/responsePseudobulk_",name_dds,"_MS_vs_CTRL.tsv"))

  # LFC shrinkage
  res_MS_shr <- lfcShrink(dds2, res = res_MS, type = "ashr")

  res_MS_shr %>% data.frame() %>% rownames_to_column("gene") %>% arrange(pvalue) %>%
    write_tsv(paste0("../../out/table/ManualClean/responsePseudobulkShr_",name_dds,"_MS_vs_CTRL.tsv"))
})

# read in the results
# plot the volcanos per cluster -------------------------------------------
folder <- "../../out/table/ManualClean/"
file <- dir(folder) %>%
  str_subset(pattern = "responsePseudobulkShr_")

df_res <- lapply(file, function(x){
  test_plot <- read_tsv(paste0(folder,x))
}) %>%
  setNames(file) %>%
  bind_rows(.id = "file") %>%
  mutate(treat = str_extract(file,pattern = "MS_vs_CTRL|NOLD_vs_control|NOLD_vs_MCI"),
         cell_id = str_extract(file,pattern = "clust_0|clust_2|clust_5|clust_6|clust_13"))

# show the distribution of the pvalues
# df_res %>%
#   mutate(comparison = str_extract(annotation,pattern = "NOLD_vs_control|MCI_vs_control|NOLD_vs_MCI")) %>%
#   ggplot(aes(x = p_val)) + geom_histogram()+theme_bw()+facet_grid(comparison~cluster)+theme(strip.background = element_blank())
# ggsave("out/image/dist_p_value_harmony_5K_scaleSbatch_SCtypeAnnotation.pdf",width = 20,height = 10)

# render all of them as a volcano plot
test_significant <- df_res %>%
  mutate(threshold = case_when(abs(log2FoldChange) > 1 & padj<0.05~1,
                               T~0)) %>%
  filter(threshold == 1)

# library(ggrepel)
df_res %>%
  # mutate(comparison = str_extract(annotation,pattern = "NOLD_vs_control|MCI_vs_control|NOLD_vs_MCI")) %>%
  # filter(symbol %in% setdiff(GOI_SLC,GOI)) %>%
  ggplot(aes(x = log2FoldChange,y = -log(padj))) +
  geom_point(alpha = 0.01) +
  geom_point(data = test_significant,aes(x = log2FoldChange,y = -log(padj)),col="red",alpha = 0.5) +
  geom_vline(xintercept = c(-1,1),col="gray",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="gray",linetype="dashed")+
  geom_text_repel(data = test_significant,aes(x = log2FoldChange,y = -log(padj),label = gene)) +
  facet_grid(treat~cell_id) +
  theme_bw()+theme(strip.background = element_blank())
ggsave("../../out/image/ManualClean/volcano_pseudobulk_MS_vs_CTRL.pdf",width = 20,height = 10)
