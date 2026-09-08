# AIM ---------------------------------------------------------------------
# follow up on scr/40_plotGene_integrationSkip_manualClean_harmony_update.R same TSPO pos vs neg DGE approach (MG/IMM and VAS cells), but:
# 1) skip the pathway/enrichR analysis
# 2) further split each cell subset by disease status (CTRL/MS) before
#    running the pos vs neg DGE, i.e. compare pos vs neg within CTRL and
#    within MS separately, for both MG and VAS
# 3) add an upset plot (ComplexUpset) to show which DE genes are shared vs
#    unique between the CTRL and MS comparisons, per cell type

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(ComplexUpset)
library(UpSetR)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")

# confirm the dataset identity
DimPlot(data.combined,group.by = "expertAnno.l1",raster = T)

GOI <- "TSPO"

table(data.combined$disease)

DefaultAssay(data.combined) <- "RNA"

# helper: subset a cell type, tag TSPO pos/neg -----------------------------
add_TSPO_cat <- function(sobj){
  exp_TSPO <- FetchData(sobj,GOI,slot = "data") %>%
    rownames_to_column("barcodes")

  sobj$TSPO <- sobj@meta.data %>%
    rownames_to_column("barcodes") %>%
    left_join(exp_TSPO,by = "barcodes") %>%
    pull(GOI)

  sobj$TSPO_cat <- case_when(sobj$TSPO == 0~"neg",
                             T~"pos")
  sobj
}

# subset the cell types of interest and tag TSPO pos/neg
sobj_MG <- subset(data.combined,subset = expertAnno.l1 == 'IMM') %>% add_TSPO_cat()
sobj_VAS <- subset(data.combined,subset = expertAnno.l1 == 'VAS') %>% add_TSPO_cat()

list_sobj <- list(MG = sobj_MG, VAS = sobj_VAS)

# check how many cells are there per condition before running the DGE
table(sobj_MG$TSPO_cat,sobj_MG$disease)
table(sobj_VAS$TSPO_cat,sobj_VAS$disease)

# TSPO DGE pos vs neg split by disease (CTRL/MS) ---------------------------
# same DGE call as in 40_plotGene_integrationSkip_manualClean_harmony_update.R
# (RunPresto, ident.1 = "pos", ident.2 = "neg"), stratified within each
# disease group
run_DGE_TSPO_disease <- function(sobj,cell_id){
  map(sort(unique(sobj$disease)),function(disease_id){
    print(paste(cell_id,disease_id))
    sobj_sub <- subset(sobj,subset = disease == disease_id)
    Idents(sobj_sub) <- "TSPO_cat"

    table(Idents(sobj_sub))

    # avg_log2FC: positive values indicate higher expression in the pos group
    RunPresto(object = sobj_sub,ident.1 = "pos",ident.2 = "neg") %>%
      rownames_to_column("gene") %>%
      mutate(cell_id = cell_id,disease = disease_id)
  }) %>%
    bind_rows()
}

df_DGE_TSPO_disease <- imap(list_sobj,run_DGE_TSPO_disease) %>%
  bind_rows()

# save the table of DGE results
df_DGE_TSPO_disease %>%
  write_tsv("../../out/table/custom/08_TSPO_res_test_posVSneg_byDisease.tsv")

# annotate DE genes and plot a volcano per cell type / disease -------------
# same thresholds as in script 40
df_volcano_TSPO_disease <- df_DGE_TSPO_disease %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no"))

ggplot() +
  geom_point(data = df_volcano_TSPO_disease %>% filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=0.5,alpha=0.2) +
  geom_point(data = df_volcano_TSPO_disease %>% filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=0.5,alpha=0.2,col="red") +
  ggrepel::geom_text_repel(data = df_volcano_TSPO_disease %>% filter(DE_cat!="no") %>% group_by(cell_id,disease) %>% slice_min(order_by = p_val_adj,n = 15),
                           aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene),size = 2,max.overlaps = 15) +
  facet_wrap(cell_id~disease,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/custom/08_TSPO_volcano_test_posVSneg_byDisease.pdf",width = 12,height = 10)

# upset plots: genes in common between CTRL and MS, per cell type ---------
# for each cell type, build separate up- and down-regulated gene sets per disease group and compare their overlap with ComplexUpset
plot_upset_by_disease <- function(cell_id){
  df_cell <- df_volcano_TSPO_disease %>%
    filter(cell_id == !!cell_id)

  list_sig_up <- lapply(sort(unique(df_cell$disease)),function(disease_id){
    df_cell %>%
      filter(disease == disease_id,DE_cat == "up") %>%
      pull(gene)
  }) %>%
    setNames(sort(unique(df_cell$disease)))

  list_sig_down <- lapply(sort(unique(df_cell$disease)),function(disease_id){
    df_cell %>%
      filter(disease == disease_id,DE_cat == "down") %>%
      pull(gene)
  }) %>%
    setNames(sort(unique(df_cell$disease)))

  (ComplexUpset::upset(fromList(list_sig_up),colnames(fromList(list_sig_up)),wrap = T) + ggtitle(paste0(cell_id," - genes up"))) +
    (ComplexUpset::upset(fromList(list_sig_down),colnames(fromList(list_sig_down)),wrap = T) + ggtitle(paste0(cell_id," - genes down")))
}

plot_upset_by_disease2 <- function(cell_id){
  df_cell <- df_volcano_TSPO_disease %>%
    filter(cell_id == !!cell_id)
  
  list_sig_up <- lapply(sort(unique(df_cell$disease)),function(disease_id){
    df_cell %>%
      filter(disease == disease_id,DE_cat == "up") %>%
      pull(gene)
  }) %>%
    setNames(sort(unique(df_cell$disease)))
  
  list_sig_down <- lapply(sort(unique(df_cell$disease)),function(disease_id){
    df_cell %>%
      filter(disease == disease_id,DE_cat == "down") %>%
      pull(gene)
  }) %>%
    setNames(sort(unique(df_cell$disease)))
  
  (ComplexUpset::upset(fromList(list_sig_up),colnames(fromList(list_sig_up)),wrap = T) + ggtitle(paste0(cell_id," - genes up")))
}

# ComplexUpset (wrap = T) returns a patchwork object whose print method draws
# the grob directly without registering itself via ggplot2::last_plot(), so
# ggsave() must be given the plot explicitly rather than relying on the
# implicit last-plot lookup
p_upset_MG <- plot_upset_by_disease("MG")
p_upset_MG
ggsave(plot = p_upset_MG,filename = "../../out/image/custom/08_TSPO_upset_MG_posVSneg_byDisease.pdf",width = 10,height = 5)

p_upset_VAS <- plot_upset_by_disease2("VAS")
p_upset_VAS
ggsave(plot = p_upset_VAS,filename = "../../out/image/custom/08_TSPO_upset_VAS_posVSneg_byDisease.pdf",width = 5,height = 5)

# save the gene identity of each upset intersection ------------------------
# for a named list of gene sets, tag every gene with the sorted, pipe-joined names of every set it belongs to (e.g. "CTRL|MS" = gene is DE in both)
get_gene_intersections <- function(list_genes){
  # long table: one row per gene per set it belongs to
  df1 <- lapply(list_genes,function(x){
    data.frame(gene = x)
  }) %>%
    bind_rows(.id = "set")

  # every unique gene across all sets
  df2 <- data.frame(gene = unique(unlist(list_genes)))

  # for each gene, pull the (sorted) combination of sets it is found in
  lapply(df2$gene,function(x){
    intersection <- df1 %>%
      dplyr::filter(gene == x) %>%
      arrange(set) %>%
      pull(set) %>%
      paste0(collapse = "|")

    data.frame(gene = x,intersection = intersection)
  }) %>%
    bind_rows()
}

# rebuild the same up/down gene sets used by the upset plots above, per cell type
get_sig_lists_by_direction <- function(cell_id,de_direction){
  df_cell <- df_volcano_TSPO_disease %>%
    filter(cell_id == !!cell_id)

  lapply(sort(unique(df_cell$disease)),function(disease_id){
    df_cell %>%
      filter(disease == disease_id,DE_cat == de_direction) %>%
      pull(gene)
  }) %>%
    setNames(sort(unique(df_cell$disease)))
}

df_upset_intersections <- imap(list_sobj,function(sobj,cell_id){
  bind_rows(
    get_gene_intersections(get_sig_lists_by_direction(cell_id,"up")) %>% mutate(direction = "up"),
    get_gene_intersections(get_sig_lists_by_direction(cell_id,"down")) %>% mutate(direction = "down")
  ) %>%
    mutate(cell_id = cell_id)
}) %>%
  bind_rows()

# save the per-gene intersection identity (which disease group(s) each DE gene falls into)
df_upset_intersections %>%
  write_tsv("../../out/table/custom/08_TSPO_upset_intersection_genes_byDisease.tsv")

# summarise how many genes fall in each intersection, per cell type / direction
df_upset_intersections %>%
  group_by(cell_id,direction,intersection) %>%
  summarise(n = n(),.groups = "drop") %>%
  arrange(cell_id,direction,desc(n))
