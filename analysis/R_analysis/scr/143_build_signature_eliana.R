# AIM ---------------------------------------------------------------------
# pull the signature files from Eliana's analysis

# libraries ---------------------------------------------------------------
library(tidyverse)
library(GSVA)
library(limma)
library(ComplexHeatmap)
library(AnnotationHub)
library(AnnotationDbi)
library(msigdbr)
library(ggrepel)
library(biomaRt)

# load the files from the DE ----------------------------------------------
# load the data an fill the missing information
df_ms_vs_ctr <- read_csv("../../data/signatures/eliana_endo_ms/Eliana_endo_deg_ms_ctrl.csv") %>%
  mutate(id1 = "ECFC_MSvsCTRL") %>%
  mutate(direction = case_when(logFC > 0 ~"_POS",
                               T ~ "_NEG")) %>%
  mutate(id = paste0(id1,direction))

df_naive_vs_trt <- read_csv("../../data/signatures/eliana_endo_ms/Eliana_endo_deg_naive_trt.csv") %>%
  mutate(id1 = "ECFC_NAIVEvsTRT") %>%
  mutate(direction = case_when(logFC > 0 ~"_POS",
                               T ~ "_NEG")) %>%
  mutate(id = paste0(id1,direction))

# build siganture file ----------------------------------------------------

df_sig_custom <- bind_rows(list(df_ms_vs_ctr,df_naive_vs_trt)) %>%
  dplyr::rename(Pathway = id, Genes = gene) %>%
  arrange(Pathway)

df_sig_custom %>%
  split(f = .$Pathway) %>%
  saveRDS("../../out/object/143_custom_signatures_eliana_endo.rds")

# -------------------------------------------------------------------------
# explore the signatures, is there any communality between the signatures?
# make a list of the significnat genes. make sure to devide coherent up and coherent down in both dataset
list_sig <- df_sig_custom %>%
  split(f = .$Pathway) %>%
  lapply(function(x){
    x %>%
      pull(Genes)
  })

# plot option 1
# UpSetR::upset(fromList(list_sig_up), order.by = "freq") 

# plot option 2
p1 <- ComplexUpset::upset(fromList(list_sig),colnames(fromList(list_sig)),wrap=T) + ggtitle("DEGs Eliana")
ggsave(plot = p1, "../../out/image/143_upset_signature_eliana.pdf",width = 6,height = 5)

# extact the intersection
df1 <- lapply(list_sig,function(x){
    data.frame(gene = x)
  }) %>% 
    bind_rows(.id = "path")
  
# pull the inique features
df2 <- data.frame(gene=unique(unlist(list_sig)))
  
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

# save the list of the intersections
df_int %>%
  arrange(int) %>%
  write_tsv("../../out/table/143_upset_signature_eliana.tsv")

# export the list
df_int %>%
  group_by(int) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# pull the genes of the intersection
GOI <- df_int %>%
  filter(int == "ECFC_MSvsCTRL_POS|ECFC_NAIVEvsTRT_NEG")  %>%
  pull(gene)

inner_join(df_ms_vs_ctr,df_naive_vs_trt,by="gene",suffix = c(".msvsctrl",".naivevstrt")) %>%
  filter(gene %in% GOI) %>%
  ggplot(aes(x=logFC.msvsctrl,y=logFC.naivevstrt)) + geom_point() + ggrepel::geom_text_repel(aes(label = gene)) + theme_bw()
ggsave("../../out/image/143_scatter_signature_eliana.pdf",width = 8,height = 8)

df_ms_vs_ctr %>%
  filter(gene %in% GOI) %>%
  print(n = 45)

df_naive_vs_trt %>%
  filter(gene %in% GOI) %>%
  print(n = 45)
