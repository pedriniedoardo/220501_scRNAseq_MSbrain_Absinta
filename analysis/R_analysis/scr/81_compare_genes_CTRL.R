# libraries ---------------------------------------------------------------
library(tidyverse)
library(UpSetR)

# load the dataset --------------------------------------------------------
# the idea is to compare the sets of genes across the different comparisons
df_res_shr <- list(OLIGO = read_tsv("../../out/table/DE_pseudobulk_OLIGO_ctrl_refCX_shr.tsv"),
     OPC = read_tsv("../../out/table/DE_pseudobulk_OPC_ctrl_refCX_shr.tsv"),
     ASTRO = read_tsv("../../out/table/DE_pseudobulk_ASTRO_ctrl_refCX_shr.tsv"),
     MG = read_tsv("../../out/table/DE_pseudobulk_MG_ctrl_refCX_shr.tsv")) %>%
  bind_rows()

# upset plot --------------------------------------------------------------
# build a list of common elements belonging to each set of fegs
list_DE_up <- df_res_shr %>%
  split(f = .$conditionVsCX) %>%
  lapply(function(x){
    x %>%
      filter(padj < 0.05,log2FoldChange>1) %>%
      pull(symbol) %>%
      unique()
  })

glimpse(list_DE_up)

list_DE_down <- df_res_shr %>%
  split(f = .$conditionVsCX) %>%
  lapply(function(x){
    x %>%
      filter(padj < 0.05,log2FoldChange<(-1)) %>%
      pull(symbol) %>%
      unique()
  })
glimpse(list_DE_down)

# try the upset plot version
# library(UpSetR)
pdf("../../out/image/upset_DEG_UP_plot_pseudobulk_WM_vs_CX.pdf",width = 14,height = 7)
upset(fromList(list_DE_up), order.by = "freq",nsets = 7)
dev.off()

pdf("../../out/image/upset_DEG_DOWN_plot_pseudobulk_WM_vs_CX.pdf",width = 14,height = 7)
upset(fromList(list_DE_down), order.by = "freq",nsets = 7)
dev.off()

# pull the intersections
df1_UP <- lapply(list_DE_up,function(x){
  data.frame(gene = x)
}) %>%
  bind_rows(.id = "path")

df1_DOWN <- lapply(list_DE_down,function(x){
  data.frame(gene = x)
}) %>%
  bind_rows(.id = "path")

head(df1_UP)
head(df1_DOWN)

df2_UP <- data.frame(gene=unique(unlist(list_DE_up)))
df2_DOWN <- data.frame(gene=unique(unlist(list_DE_down)))

head(df2_UP)
head(df2_DOWN)

df_int_UP <- lapply(df2_UP$gene,function(x){
  # pull the name of the intersections
  intersection <- df1_UP %>%
    dplyr::filter(gene==x) %>%
    arrange(path) %>%
    pull("path") %>%
    paste0(collapse = "|")

  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>%
  bind_rows()

df_int_DOWN <- lapply(df2_DOWN$gene,function(x){
  # pull the name of the intersections
  intersection <- df1_DOWN %>%
    dplyr::filter(gene==x) %>%
    arrange(path) %>%
    pull("path") %>%
    paste0(collapse = "|")

  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>%
  bind_rows()

df_int_UP %>%
  write_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_WM_vs_CX.tsv")

df_int_DOWN %>%
  write_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_WM_vs_CX.tsv")

head(df_int_UP,n=20)
head(df_int_DOWN,n=20)

df_int_UP %>%
  group_by(int) %>%
  summarise(n=n()) %>%
  arrange(desc(n))

df_int_DOWN %>%
  group_by(int) %>%
  summarise(n=n()) %>%
  arrange(desc(n))

# explore the genes -------------------------------------------------------
# explore the meaning of the genes in common
# try to run a pathway analysis on the common genes
df_int_UP <- read_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_WM_vs_CX.tsv")
df_int_DOWN <- read_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_WM_vs_CX.tsv")

# for the genes up in WM focus on the cell type specific genes alone
list_UP_WM <- df_int_UP %>%
  dplyr::filter(int %in% c("ASTRO.WM.Ctrl_shr","MG.WM.Ctrl_shr","OPC.WM.Ctrl_shr","OLIGO.WM.Ctrl_shr")) %>%
  mutate(test = case_when(int == "ASTRO.WM.Ctrl_shr"~"ASTRO_spec_WM",
                          int == "MG.WM.Ctrl_shr"~"MG_spec_WM",
                          int == "OPC.WM.Ctrl_shr"~"OPC_spec_WM",
                          int == "OLIGO.WM.Ctrl_shr"~"OLIGO_spec_WM")) %>%
  split(f=.$test)

# for the genes up in CX select the sepcific set of genes and the one common across all
list_UP_CX <- df_int_DOWN %>%
  dplyr::filter(int %in% c("ASTRO.WM.Ctrl_shr","MG.WM.Ctrl_shr","OPC.WM.Ctrl_shr","OLIGO.WM.Ctrl_shr","ASTRO.WM.Ctrl_shr|MG.WM.Ctrl_shr|OLIGO.WM.Ctrl_shr|OPC.WM.Ctrl_shr")) %>%
  mutate(test = case_when(int == "ASTRO.WM.Ctrl_shr"~"ASTRO_spec_CX",
                          int == "MG.WM.Ctrl_shr"~"MG_spec_CX",
                          int == "OPC.WM.Ctrl_shr"~"OPC_spec_CX",
                          int == "OLIGO.WM.Ctrl_shr"~"OLIGO_spec_CX",
                          int == "ASTRO.WM.Ctrl_shr|MG.WM.Ctrl_shr|OLIGO.WM.Ctrl_shr|OPC.WM.Ctrl_shr"~"ALL_CX")) %>%
  split(f=.$test)

# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
# filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Azimuth"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Tabula"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","HDSigDB_Human_2021","GO_Biological_Process_2023","Azimuth_2023","CellMarker_2024","Tabula_Sapiens")

# query -------------------------------------------------------------------
# seelct only the clusters with more than 10 genes as degs
# file <- dir("../../out/table/") %>%
#   str_subset(pattern = "res_") %>%
#   str_subset(pattern = ".txt") %>%
#   str_subset(pattern = "HUVEC") %>%
#   str_subset(pattern = "shr",negate = F)
# file

# load the results
# results <- lapply(paste0("../../out/table/",file),function(x){
#   read_tsv(x)  %>%
#     filter(padj<0.05 & abs(log2FoldChange)>1&!is.na(symbol))
# }) %>%
#   setNames(str_remove_all(file,pattern = ".txt"))
list_genes_UP <- list_UP_WM %>%
  map(function(x){
    x %>%
      pull(gene)
  })

list_genes_DOWN <- list_UP_CX %>%
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
  write_tsv("../../out/table/enrichR_DE_intersection_WM.tsv")

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
  write_tsv("../../out/table/enrichR_DE_intersection_CX.tsv")

# see the top term from azimuth cell type
list_enrichr_UP %>%
  as.tibble() %>%
  # group_by(annotation) %>% summarise()
  filter(annotation == "MSigDB_Hallmark_2020",
         comparison == "ASTRO_spec_WM") %>%
  # select(Term,Overlap,Genes,Adjusted.P.value,Odds.Ratio) %>%
  data.frame()  %>%
  head()

list_enrichr_UP %>%
  as.tibble() %>%
  # group_by(annotation) %>% summarise()
  filter(annotation == "MSigDB_Hallmark_2020",
         comparison == "OLIGO_spec_WM") %>%
  # select(Term,Overlap,Genes,Adjusted.P.value,Odds.Ratio) %>%
  data.frame()  %>%
  head()

list_enrichr_DOWN %>%
  as.tibble() %>%
  # group_by(annotation) %>% summarise()
  filter(annotation == "HDSigDB_Human_2021",
         comparison == "ALL_CX") %>%
  # select(Term,Overlap,Genes,Adjusted.P.value,Odds.Ratio) %>%
  data.frame()  %>%
  head()

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
ggsave("../../out/image/enrichR_DE_intersection_WM.pdf",width = 24,height = 20,limitsize = FALSE)

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
ggsave("../../out/image/enrichR_DE_intersection_CX.pdf",width = 30,height = 20,limitsize = FALSE)

# -------------------------------------------------------------------------
read_tsv("../../out/table/DE_pseudobulk_OLIGO_ctrl_refCX_shr.tsv") |>
  filter(symbol %in% c("HEXB"))
