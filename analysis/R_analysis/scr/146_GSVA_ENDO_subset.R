# AIM ---------------------------------------------------------------------
# run GSVA on the ENDO subset of the brain daataset

# libraries ---------------------------------------------------------------
library(tidyverse)
library(GSVA)
library(limma)
library(ComplexHeatmap)
# library(org.Hs.eg.db)
library(AnnotationHub)
library(AnnotationDbi)
library(msigdbr)
library(ggrepel)

# read in the data --------------------------------------------------------
# read in the filtered pbulk object for MG subset
ddsHTSeq_filter <- readRDS("../../out/object/145_ddsHTSeq_pseudobulk_filterExp_ENDO.rds")

exp <- counts(ddsHTSeq_filter,normalized = T) %>%
  data.frame() %>%
  # dplyr::select(contains("BASELINE")|contains("Fe")|contains("myelin")) %>%
  as.matrix()

# read in the signature file provided by Eliana
gene_sets <- readRDS("../../out/object/143_custom_signatures_eliana_endo.rds")

# format in order to be accepted by GSEA
pathways <- lapply(gene_sets,function(x){
  x %>%
    pull(Genes)
})

# -------------------------------------------------------------------------
# set.seed for reproducibility
set.seed(123)

# perform the GSEA based on the normalized counts
es <- gsva(exp,
           pathways,
           min.sz=2,
           max.sz=500,
           kcdf="Poisson",
           mx.diff=TRUE,
           verbose=FALSE,
           parallel.sz=1)

es_log <- gsva(log(exp+1),
               pathways,
               min.sz=2,
               max.sz=500,
               kcdf="Gaussian",
               mx.diff=TRUE,
               verbose=FALSE,
               parallel.sz=1)

# show correlation between the two estimates ------------------------------
left_join(
  es %>%
    data.frame() %>%
    rownames_to_column("pathway") %>%
    pivot_longer(names_to = "sample",values_to = "estimate",-pathway),
  es_log %>%
    data.frame() %>%
    rownames_to_column("pathway") %>%
    pivot_longer(names_to = "sample",values_to = "estimate_log",-pathway),by = c("pathway","sample")) %>%
  ggplot(aes(x=estimate_log,y=estimate))+geom_point()

# STATISTICAL TESTING -----------------------------------------------------
# run the analysis only on the logged normalized values

# use limma for testing significance
# library(limma)
# adjPvalueCutoff <- 0.001
# logFCcutoff <- log2(2)
# logFCcutoff

# define the factor for the treatmentnt basesd on the colnames
lut <- data.frame(colname = colnames(es)) %>%
  mutate(origin = unlist(str_extract_all(colname,pattern = "cortex|wm"))) %>%
  mutate(disease = unlist(str_extract_all(colname,pattern = "CTRL|MS")))

design <- model.matrix(~ lut$origin + lut$disease)
colnames(design)[1] <- c("intercept")
design

# fit <- lmFit(es, design)
fit_log <- lmFit(es_log, design)
# fit <- eBayes(fit)
fit_log <- eBayes(fit_log)


# allGenesets_Fe <- topTable(fit, coef="lut$treatFe", number=Inf)
allGenesets_log <- topTable(fit_log, coef="lut$diseaseMS", number=Inf)

# pull the stats from Aletta's list
allGenesets_log %>%
  rownames_to_column("GO_id")

# save the tables
allGenesets_log %>%
  rownames_to_column("GO_id") %>%
  left_join(test01_01_summary,by = c("GO_id"="gs_exact_source"))
  # write_tsv("../../out/table/24_df_table_GSVA_GOBP_Fe_MG.tsv")

# volcano GSVA ------------------------------------------------------------
df_plot <- allGenesets_log %>%
  rownames_to_column("pathway")

df_plot2 <- df_plot %>%   
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "HALLMARK_|KEGG_|GOBP_") %>%
           str_sub(start = 1,end = 35))
# mutate(color = case_when(pathway %in% c("KEGG_PRIMARY_IMMUNODEFICIENCY",
#                                       "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
#                                       "KEGG_VEGF_SIGNALING_PATHWAY",
#                                       "KEGG_CELL_CYCLE",
#                                       "KEGG_CELL_ADHESION_MOLECULES_CAMS",
#                                       "KEGG_ECM_RECEPTOR_INTERACTION")~"red",
#                        T~"black")) %>%
# mutate(color = factor(color))

df_plot2 %>%
  # ggplot(aes(y = -log10(adj.P.Val),x = logFC,label = pathway2,col=color)) +
  ggplot(aes(y = -log10(adj.P.Val),x = logFC)) +
  geom_point(alpha = 0.2) +
  # geom_point(aes(size = size),alpha = 0.2) +
  # facet_wrap(~dataset) +
  # theme_bw(base_rect_size = 2)+
  theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA),
  #       axis.ticks = element_line(colour = "black"),
  #       #axis.ticks.length = unit(.25, "cm")
  #       legend.position = "none"
  # )+
  theme(strip.background = element_blank(),legend.position = "none") +
  geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray",alpha=0.8)

# plot heatmap ------------------------------------------------------------
# library(ComplexHeatmap)
# define only the significnat terms
DEgeneSets <- allGenesets_log_Fe %>%
  rownames_to_column("pathway") %>%
  # filter(adj.P.Val < 0.05) %>%
  pull(pathway) 

# define the matrix for the heatmap
mat_norm <- es %>%
  data.frame() %>%
  rownames_to_column() %>%
  # plot only the significant terms
  filter(rowname %in% DEgeneSets) %>%
  # scale the values rowwise
  gather(key = sample,value = exp,-rowname) %>%
  group_by(rowname) %>%
  mutate(norm = (exp - mean(exp))/sd(exp)) %>%
  dplyr::select(-exp) %>%
  spread(key = sample,value = norm) %>%
  column_to_rownames()

sample_ordered <- str_extract(colnames(mat_norm),pattern = "CTRL|MS")
sample_ordered

# build the annotation object
column_ha <- HeatmapAnnotation(treat = sample_ordered,
                               col = list(treat = c("CTRL" = "blue", "MS"="red"))) 

hm <- Heatmap(mat_norm,show_row_names = F,
                 # add annotation for the columns
                 # hide columns labels
                 # show_column_names = F,
                 # fix width of the lables
                 top_annotation = column_ha,
                 # row_names_gp = gpar(col = rowname_color),
                 row_names_max_width = max_text_width(
                   rownames(mat_norm),
                   gp = gpar(fontsize = 12)
                 ))

# pdf(file = "../../out/image/24_heatmap_GSVA_GOBP_Fe_es_Log_ZScore.pdf", width = 8, height = 5)
draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 100), "mm"))
# dev.off()
