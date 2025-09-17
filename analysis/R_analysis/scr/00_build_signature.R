# AIM ---------------------------------------------------------------------
# build the signature object for the targeted analysis by GSEA

# libraries ---------------------------------------------------------------
library(tidyverse)
library(GSEABase)

# read in the gmt siganture files -----------------------------------------
# read in the autophagy set of signatures
file_gmt <- dir("../../data/signatures/autophagy/source/") %>%
  str_subset(pattern = ".gmt")

name_gmt <- file_gmt %>%
  str_sub(start = 1,end = -16)

pathways_gmt <- lapply(file_gmt, function(x){
  name <- paste0("../../data/signatures/autophagy/source/",x)
  getGmt(name) %>%
    geneIds() %>%
    .[[1]]
}) %>%
  setNames(name_gmt)

# rearrange the shape of the signature
pathways_gmt2 <- pmap(list(pathways_gmt,names(pathways_gmt)),function(x,name){
  df <- data.frame(Pathway = name,Genes = x)
  return(df)
})

# # read in the LUT to convert the gene names form mouse to human
# homologeneData2_240528 <- read_tsv("../../data/Homologene_240528") %>%
#   as.data.frame()
# 
# s.genes.mouse <- homologene(s.genes, inTax = 9606, outTax = 10090,db = homologeneData2_240528) %>%
#   pull(`10090`) %>%
#   unique()
# g2m.genes.mouse <- homologene(g2m.genes, inTax = 9606, outTax = 10090,db = homologeneData2_240528) %>%
#   pull(`10090`) %>%
#   unique()

# # read in the txt signature files -----------------------------------------
# # the review sigantures
# file_txt <- dir("../../data/signatures/signatures/") %>%
#   str_subset(pattern = "review.txt")
# 
# name_txt <- file_txt %>%
#   str_sub(start = 1,end = -12)
# 
# pathways_txt <- lapply(file_txt, function(x){
#   name <- paste0("../../data/signatures/",x)
#   read_tsv(name) %>% 
#     dplyr::select(human_gene) %>% 
#     drop_na() %>%
#     pull(human_gene)
#   
# }) %>%
#   setNames(name_txt)
# 
# # read in the txt signatures file (not review) ----------------------------
# # the non review signatures
# file_txt2 <- dir("../../data/signatures/") %>%
#   str_subset(pattern = ".txt") %>%
#   str_subset(pattern = "review",negate = T)
# 
# name_txt2 <- file_txt2 %>%
#   str_sub(start = 1,end = -5)
# 
# pathways_txt2 <- lapply(file_txt2, function(x){
#   name <- paste0("../../data/signatures/",x)
#   read_tsv(name) %>% 
#     dplyr::select(human_gene) %>% 
#     drop_na() %>%
#     pull(human_gene)
#   
# }) %>%
#   setNames(name_txt2)
# 
# # read in the csv signatures file (not review) ----------------------------
# # the non review signatures
# 
# pathways_csv_all <- read_csv("../../data/signatures/marker_vs_all.csv") %>%
#   mutate(Cell_Type = str_remove_all(Cell_Type,pattern = "BEC, ")) %>%
#   dplyr::filter(Cell_Type %in% c("Arterial","Capillary","Venous")) %>%
#   mutate(Cell_Type = paste0(Cell_Type,"_vs_all")) %>%
#   split(f = .$Cell_Type) %>%
#   map(function(x){
#     x %>% pull(Gene)
#   })
# 
# pathways_csv_mural <- read_csv("../../data/signatures/marker_vs_mural.csv") %>%
#   # group_by(Cell_subtype) %>%
#   # summarise()
#   mutate(Cell_subtype = str_remove_all(Cell_subtype,pattern = "-like/ Proteostatic")) %>%
#   mutate(Cell_subtype = paste0(Cell_subtype,"_vs_mural")) %>%
#   split(f = .$Cell_subtype) %>%
#   map(function(x){
#     x %>% pull(Gene)
#   })


# -------------------------------------------------------------------------
# merge all the pathways in a single list
# pathways <- c(pathways_txt2,
#               pathways_txt,
#               pathways_gmt,
#               pathways_csv_mural,
#               pathways_csv_all)

pathways <- c(pathways_gmt2)

# save the object
saveRDS(pathways,"../../data/signatures/autophagy_pathways.rds")
