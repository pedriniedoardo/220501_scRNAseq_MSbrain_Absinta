# libraries ---------------------------------------------------------------
library(tidyverse)
library(UpSetR)
library(cowplot)
library(enrichR)
library(scales)
library(patchwork)

# load the dataset --------------------------------------------------------
# the idea is to compare the sets of genes across the different comparisons

# list of comparisons for MS vs CTRL in WM
df_res_shr_WM <- list(OLIGO = read_tsv("../../out/table/DE_pseudobulk_OLIGO_full_WM_shr.tsv"),
                      OPC = read_tsv("../../out/table/DE_pseudobulk_OPC_full_WM_shr.tsv"),
                      ASTRO = read_tsv("../../out/table/DE_pseudobulk_ASTRO_full_WM_shr.tsv"),
                      MG = read_tsv("../../out/table/DE_pseudobulk_MG_WM_shr.tsv")) %>%
  bind_rows(.id = "cell_id") %>%
  mutate(tissue = "WM")

# list of comparisons for MS vs CTRL in CX
df_res_shr_CX <- list(OLIGO = read_tsv("../../out/table/DE_pseudobulk_OLIGO_full_CX_shr.tsv"),
                      OPC = read_tsv("../../out/table/DE_pseudobulk_OPC_full_CX_shr.tsv"),
                      ASTRO = read_tsv("../../out/table/DE_pseudobulk_ASTRO_full_CX_shr.tsv"),
                      MG = read_tsv("../../out/table/DE_pseudobulk_MG_CX_shr.tsv")) %>%
  bind_rows(.id = "cell_id") %>%
  mutate(tissue = "CX")

# plot gene count ---------------------------------------------------------
# plot the counts of genes crossing the clessical threshold of significance 
bind_rows(df_res_shr_CX,df_res_shr_WM) %>%
  mutate(sig = case_when(padj<0.05 & abs(log2FoldChange)>1 ~ "sig",
                         T~"non_sig")) %>%
  mutate(direction = sign(log2FoldChange)) %>%
  mutate(direction = factor(direction)) %>%
  filter(sig == "sig") %>%
  group_by(cell_id,tissue,direction) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(cell_id) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(cell_id = fct_reorder(cell_id,tot,.desc = T)) %>%
  ggplot(aes(x=tissue,y=n,fill=direction))+
  geom_col()+
  facet_wrap(~cell_id,nrow = 1)+
  theme_cowplot()+
  theme(strip.background = element_blank())+scale_fill_manual(values = c("blue","red"))
ggsave("../../out/image/barplot_DEG_MS_vs_CTRL_all.pdf",width = 6,height = 3)

# plot volcano fro all comparisons ----------------------------------------
plot_volcano <- bind_rows(df_res_shr_CX,df_res_shr_WM) %>%
  mutate(sig = case_when(padj<0.05 & abs(log2FoldChange)>1 ~ "sig",
                         T~"non_sig")) %>%
  mutate(direction = sign(log2FoldChange)) %>%
  mutate(direction = factor(direction)) %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # ggrepel::geom_text_repel( data = plot_volcano[plot_volcano$col==1,] %>% group_by(conditionVsCTRL) %>% arrange(padj) %>% dplyr::slice(1:10), aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(tissue~cell_id,ncol=4,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/image/volcano_DEG_MS_vs_CTRL_all.pdf",width = 12,height = 6)

# upset plot --------------------------------------------------------------
# considering the situation of the dataset I will be working only with WM data
# build a list of common elements belonging to each set of fegs
list_DE_up <- df_res_shr_WM %>%
  split(f = .$conditionVsCTRL) %>%
  lapply(function(x){
    x %>%
      filter(padj < 0.05,log2FoldChange>1) %>%
      pull(symbol) %>%
      unique()
  })

glimpse(list_DE_up)

list_DE_down <- df_res_shr_WM %>%
  split(f = .$conditionVsCTRL) %>%
  lapply(function(x){
    x %>%
      filter(padj < 0.05,log2FoldChange<(-1)) %>%
      pull(symbol) %>%
      unique()
  })
glimpse(list_DE_down)

# try the upset plot version
# library(UpSetR)
pdf("../../out/image/upset_DEG_UP_plot_pseudobulk_MS_vs_CTRL_WM.pdf",width = 14,height = 7)
upset(fromList(list_DE_up), order.by = "freq",nsets = 7)
dev.off()

pdf("../../out/image/upset_DEG_DOWN_plot_pseudobulk_MS_vs_CTRL_WM.pdf",width = 14,height = 7)
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
  bind_rows() %>%
  arrange(int)

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
  bind_rows() %>%
  arrange(int)

df_int_UP %>%
  write_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_MS_vs_CTRL_WM.tsv")

df_int_DOWN %>%
  write_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_MS_vs_CTRL_WM.tsv")

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

# what are the common 21 genes
df_int_UP %>%
  filter(int == "ASTRO.WM.MS_shr|MG.WM.MS_shr|OLIGO.WM.MS_shr|OPC.WM.MS_shr")

# explore the genes -------------------------------------------------------
# explore the meaning of the genes in common
# try to run a pathway analysis on the common genes
df_int_UP <- read_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_MS_vs_CTRL_WM.tsv")
df_int_DOWN <- read_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_MS_vs_CTRL_WM.tsv")

# for the genes up in WM focus on the cell type specific genes alone
list_UP_WM <- df_int_UP %>%
  dplyr::filter(int %in% c("ASTRO.WM.MS_shr|MG.WM.MS_shr|OLIGO.WM.MS_shr|OPC.WM.MS_shr","ASTRO.WM.MS_shr","MG.WM.MS_shr","OPC.WM.MS_shr","OLIGO.WM.MS_shr")) %>%
  mutate(test = case_when(int == "ASTRO.WM.MS_shr"~"ASTRO_spec_WM.MS",
                          int == "MG.WM.MS_shr"~"MG_spec_WM.MS",
                          int == "OPC.WM.MS_shr"~"OPC_spec_WM.MS",
                          int == "OLIGO.WM.MS_shr"~"OLIGO_spec_WM.MS",
                          int == "ASTRO.WM.MS_shr|MG.WM.MS_shr|OLIGO.WM.MS_shr|OPC.WM.MS_shr"~"ALL_WM.MS")) %>%
  split(f=.$test)

list_DOWN_WM <- df_int_DOWN %>%
  dplyr::filter(int %in% c("ASTRO.WM.MS_shr","MG.WM.MS_shr","OPC.WM.MS_shr","OLIGO.WM.MS_shr")) %>%
  mutate(test = case_when(int == "ASTRO.WM.MS_shr"~"ASTRO_spec_WM.MS",
                          int == "MG.WM.MS_shr"~"MG_spec_WM.MS",
                          int == "OPC.WM.MS_shr"~"OPC_spec_WM.MS",
                          int == "OLIGO.WM.MS_shr"~"OLIGO_spec_WM.MS")) %>%
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
list_genes_UP <- list_UP_WM %>%
  map(function(x){
    x %>%
      pull(gene)
  })

list_genes_DOWN <- list_DOWN_WM %>%
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
  write_tsv("../../out/table/enrichR_DE_intersection_UP_MS_vs_CTRL_WM.tsv")

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
  write_tsv("../../out/table/enrichR_DE_intersection_DOWN_MS_vs_CTRL_WM.tsv")

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
ggsave("../../out/image/enrichR_DE_intersection_UP_MS_vs_CTRL_WM.pdf",width = 30,height = 20,limitsize = FALSE)

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
ggsave("../../out/image/enrichR_DE_intersection_DOWN_MS_vs_CTRL_WM.pdf",width = 24,height = 20,limitsize = FALSE)
