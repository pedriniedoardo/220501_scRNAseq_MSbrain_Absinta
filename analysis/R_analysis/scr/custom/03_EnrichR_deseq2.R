# librarues ---------------------------------------------------------------
library(tidyverse)
library(enrichR)
library(scales)

# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Azimuth"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_Pathways_2024")

# GENE SELECTION ----------------------------------------------------------
# list of genes to consider for the enrichment analysis
test_out <- read_tsv("../out/table/02_DE_pseudobulk_filterExp_shr.tsv") %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > 0.5) %>%
  split(f = .$condVsCTRL)
glimpse(test_out)

# read in the metadata for the genes
# meta <- read_csv("data/scTrem/GSE130626_gene_info.csv")

#
list_res_tot <- lapply(test_out, function(x){
  x %>%
    pull(symbol)
})
# query -------------------------------------------------------------------
list <- lapply(list_res_tot,function(x){
  genes <- x
  out_enrich <- enrichr(genes, dbs_db)
  #
  out_enrich %>%
    bind_rows(.id = "annotation")
})

df_enrichr_annotation_enriched_tot <- list %>%
  bind_rows(.id = "comparison")

df_enrichr_annotation_enriched_tot %>%
  write_tsv("../out/table/03_enrichR_deseq2.tsv")

# check entries for senescence
df_enrichr_annotation_enriched_tot %>%
  filter(str_detect(Term,pattern = "enescence"))

# library(scales)
list_plot <- lapply(unique(df_enrichr_annotation_enriched_tot$comparison),function(x){
  df_enrichr_annotation_enriched_tot %>%
    filter(comparison == x) %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:20) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    theme(strip.background = element_blank()) +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))
  
}) %>%
  setNames(unique(df_enrichr_annotation_enriched_tot$comparison))

# scale_color_gradient(low = "red",high = "blue")
pmap(list(list_plot,names(list_plot)), function(x,y){
  ggsave(plot = x,paste0("../out/plot/03_enrichR_deseq2_",y,".pdf"),width = 7,height = 10)
})

#

