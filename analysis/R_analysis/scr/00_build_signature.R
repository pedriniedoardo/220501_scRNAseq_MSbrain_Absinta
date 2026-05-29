# AIM ---------------------------------------------------------------------
# build the signature object for the targeted analysis by GSEA

# libraries ---------------------------------------------------------------
library(tidyverse)
library(GSEABase)
library(msigdbr)

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

# Sofia Sifgnatures -------------------------------------------------------
# custom signatures she pulled from the databases

# manually compile the list from the table sofia has provided from MitoCarta3.0 mitochondrial related
# read in the table and split by signature
df_MitoCarta <- read_csv("../../data/signatures/Sofia/MitoCarta3_Pathways.csv") %>%
  mutate(pathway_fix = paste0("MitoCart_",str_replace_all(string = Pathway,pattern = " ",replacement = ".")))

pathways_mitocarta <- split(x = df_MitoCarta$Genes, f = df_MitoCarta$pathway_fix)

# save the object
saveRDS(pathways_mitocarta,"../../data/signatures/Sofia/mitocarta_related_pathways.rds")

# msigdb sigantures
msigdbr_collections() %>% print(n=30)
gene_sets <- msigdbr(species = "Homo sapiens")
head(gene_sets)

# inflammatory related
# BIOCARTA_IFNA_PATHWAY
# WP_CYTOKINES_AND_INFLAMMATORY_RESPONSE
# WP_TNFALPHA_SIGNALING_PATHWAY
gene_sets_inf <- gene_sets %>%
  filter(gs_id %in% c("M39662",
                      "M22056",
                      "M39711")) 
# confirm the signatures
gene_sets_inf %>%
  pull(gs_name) %>%
  unique()

# format in order to be accepted by GSEA
pathways_inf <- split(x = gene_sets_inf$gene_symbol, f = gene_sets_inf$gs_name)

# save the object
saveRDS(pathways_inf,"../../data/signatures/Sofia/inflammation_related_pathways.rds")

# cellular iron
# "GOBP_CELLULAR_IRON_ION_HOMEOSTASIS" "GOBP_CELLULAR_RESPONSE_TO_IRON_ION"                                     
# "GOBP_IRON_IMPORT_INTO_CELL" "GOBP_IRON_ION_IMPORT_ACROSS_PLASMA_MEMBRANE"                            
# "GOBP_REGULATION_OF_IRON_ION_TRANSPORT" "GOBP_RESPONSE_TO_IRON_II_ION"                                           
# "GOBP_RESPONSE_TO_IRON_ION" "GOCC_IRON_SULFUR_CLUSTER_ASSEMBLY_COMPLEX"                              
# "GOMF_IRON_ION_BINDING" "GOMF_IRON_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY"                       
# "HP_IRON_ACCUMULATION_IN_BRAIN" "WP_NEURODEGENERATION_WITH_BRAIN_IRON_ACCUMULATION_NBIA_SUBTYPES_PATHWAY"
gene_sets_iron <- gene_sets %>%
  filter(gs_id %in% c("M24450",
                      "M11074",
                      "M40373",
                      "M16655",
                      "M43154",
                      "M18915",
                      "M23300",
                      "M34268",
                      "M46774",
                      "M23194",
                      "M26259",
                      "M38155",
                      "M39602")) 
# confirm the signatures
gene_sets_iron %>%
  pull(gs_name) %>%
  unique()

# format in order to be accepted by GSEA
pathways_iron <- split(x = gene_sets_iron$gene_symbol, f = gene_sets_iron$gs_name)

# save the object
saveRDS(pathways_iron,"../../data/signatures/Sofia/iron_related_pathways.rds")


# lipid related
# "BIOCARTA_ARENRF2_PATHWAY" "GOBP_C21_STEROID_HORMONE_BIOSYNTHETIC_PROCESS" "GOBP_GLUTATHIONE_METABOLIC_PROCESS" "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION"    
# "WP_NRF2_PATHWAY"
gene_sets_lipid <- gene_sets %>%
  filter(gs_id %in% c("M49532",
                      "M15204",
                      "M46862",
                      "M50529",
                      "M14708",
                      "M14339",
                      "M39454",
                      "M2164",
                      "M19236"
  )) 
# confirm the signatures
gene_sets_lipid %>%
  pull(gs_name) %>%
  unique()

# format in order to be accepted by GSEA
pathways_lipid <- split(x = gene_sets_lipid$gene_symbol, f = gene_sets_lipid$gs_name)

# save the object
saveRDS(pathways_lipid,"../../data/signatures/Sofia/lipid_related_pathways.rds")


# metabolism
# "GOBP_CELLULAR_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES" "GOBP_OXIDATIVE_PHOSPHORYLATION"                                                            
# "GOBP_POSITIVE_REGULATION_OF_OXIDATIVE_STRESS_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY" "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"
gene_sets_metabolism <- gene_sets %>%
  filter(gs_id %in% c("M11521",
                      "M12919",
                      "M45123",
                      "M25058",
                      "M16581")) 
# confirm the signatures
gene_sets_metabolism %>%
  pull(gs_name) %>%
  unique()

# format in order to be accepted by GSEA
pathways_metabolism <- split(x = gene_sets_metabolism$gene_symbol, f = gene_sets_metabolism$gs_name)

# save the object
saveRDS(pathways_metabolism,"../../data/signatures/Sofia/metabolism_related_pathways.rds")
