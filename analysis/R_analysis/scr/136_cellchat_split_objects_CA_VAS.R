# AIM ---------------------------------------------------------------------
# this script will split the object that martina as requested. in particular we are focussing on the CA samples only. In one object put all the cells from CA and the VAS cells called as senescent, in the other object the same cells, but the VAS not called as senescent

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(circlize)
library(ComplexHeatmap)
library(CellChat)
library(UpSetR)

# read in the data --------------------------------------------------------
# pull the full object
Seurat.object <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
DimPlot(Seurat.object,raster = T,group.by = "expertAnno.l1")

# add the barcode variable to the metadata to make it easier the filtering
Seurat.object$barcodes <- colnames(Seurat.object)
Seurat.object@meta.data

# pull the annotation from SIT
LUT_SIT <- read_tsv("../../out/table/revision/124_meta_SIT_WM_CX_harmonySkipIntegration_AllSoupX_test.tsv")

# how many cells from the IMM cluster in CA are labelled as senescent
LUT_SIT %>%
  filter(pathology_class == "WM_CA") %>%
  group_by(expertAnno.l1,SENEQUANTILE) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = SENEQUANTILE,values_from = n) %>%
  mutate(tot = NO + YES)

# pull the annotation from Senmayo method
# folder <-  "../../out/table/revision/modules_SENESCENCE/"
# # focus on the autophagy only sigantures
# pattern_sig <- paste0(names(readRDS("../../data/signatures/senescence_pathways.rds")),collapse = "|")
# 
# file <- dir(folder) %>% 
#   str_subset(pattern = "Module") %>%
#   str_subset(pattern = pattern_sig)
# 
# df_modules_cellID <- lapply(file, function(x){
#   read_tsv(paste0(folder,x))
# }) %>% 
#   bind_rows() %>% 
#   mutate(disease = factor(disease,levels = c("CTRL","MS")))
# 
# LUT_sample <- read_tsv("../../out/table/revision/120_meta.data_AnnotationSCType_manualAnnotation.tsv") %>% 
#   group_by(orig.ident,pathology_class) %>% 
#   summarise()
# 
# # use wm ctrl as reference
# id_signature <- unique(df_modules_cellID$signature)
# 
# df_90_cellID2 <- df_modules_cellID %>% 
#   filter(signature %in% "senmayo",
#          disease %in% "CTRL",
#          origin == "cortex") %>% 
#   group_by(signature,expertAnno.l1) %>%
#   summarise(thr = quantile(signature_score1,prob=0.90))
#   
# test_summary_cellID <- df_modules_cellID %>% 
#   filter(signature %in% "senmayo") %>%
#   left_join(df_90_cellID2,c("expertAnno.l1","signature")) %>% 
#   mutate(sen = case_when(signature_score1 > thr~1,
#                          T~0))
# 
# test_summary_cellID %>%
#   write_tsv("../../out/table/revision/124_meta_senmayo_WM_CX_harmonySkipIntegration_AllSoupX_test.tsv")

LUT_senmayo <- read_tsv("../../out/table/revision/124_meta_senmayo_WM_CX_harmonySkipIntegration_AllSoupX_test.tsv")

# how many cells from the IMM cluster in CA are labelled as senescent
LUT_senmayo %>%
  filter(pathology_class == "WM_CA") %>%
  group_by(expertAnno.l1,sen) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = sen,values_from = n) %>%
  mutate(tot = `0`+`1`)

# check if this cells are the same across senmayo and SIT call
# list_barcodes <- list(senmayo = LUT_senmayo %>%
#                         filter(pathology_class == "WM_CA",
#                                expertAnno.l1 == 'IMM',sen == 1) %>%
#                         pull(barcode),
#                       SIT = LUT_SIT %>%
#                         filter(pathology_class == "WM_CA",
#                                expertAnno.l1 == 'IMM',SENEQUANTILE == "YES") %>%
#                         pull(rowname))

# upset(fromList(list_barcodes), order.by = "freq")
list_barcodes <- list(senmayo = LUT_senmayo %>%
                        filter(pathology_class == "WM_CA",
                               expertAnno.l1 == 'VAS',sen == 1) %>%
                        pull(barcode),
                      SIT = LUT_SIT %>%
                        filter(pathology_class == "WM_CA",
                               expertAnno.l1 == 'VAS',SENEQUANTILE == "YES") %>%
                        pull(rowname))


# pdf("../../out/image/revision/126_upset_senescence_signature.pdf",width = 5,height = 3,onefile = T)
# upset(fromList(list_barcodes), order.by = "freq")
# dev.off()
plot <- ComplexUpset::upset(fromList(list_barcodes),colnames(fromList(list_barcodes)),wrap=T) + ggtitle("VAS senescence")
ggsave(plot = plot,"../../out/image/revision/134_upset_senescence_signature_VAS.pdf",width = 5,height = 3)

# df1 <- lapply(list_filter,function(x){
#   data.frame(gene = x)
# }) %>% 
#   bind_rows(.id = "path")
# 
# head(df1)
# 
# df2 <- data.frame(gene=unique(unlist(list_filter)))
# 
# head(df2)
# 
# df_int <- lapply(df2$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1 %>% 
#     dplyr::filter(gene==x) %>% 
#     arrange(path) %>% 
#     pull("path") %>% 
#     paste0(collapse = "|")
#   
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>% 
#   bind_rows()
# 
# head(df_int,n=20)
# 
# df_int %>% 
#   group_by(int) %>% 
#   summarise(n=n()) %>% 
#   arrange(desc(n))

# split the object --------------------------------------------------------
# senmayo
# pull all the cells from CA sample. for the IMM subset pull either the non senescent or the senescent cells
barcodes_senSenmayo <- LUT_senmayo %>%
  filter(pathology_class == "WM_CA") %>%
  mutate(keep = case_when(
    expertAnno.l1 != "VAS" ~ 1,
    expertAnno.l1 == "VAS" & sen == 1 ~ 1,
    T~0)) %>%
  filter(keep ==1) %>%
  pull(barcode)

barcodes_ctrSenmayo <- LUT_senmayo %>%
  filter(pathology_class == "WM_CA") %>%
  mutate(keep = case_when(
    expertAnno.l1 != "VAS"~1,
    expertAnno.l1 == "VAS" & sen == 0~1,
    T~0)) %>%
  filter(keep ==1) %>%
  pull(barcode)

length(barcodes_ctrSenmayo)
length(barcodes_senSenmayo)

# 3019+182+1010+2864+92+217+6541+570+
# 632
# 159

# save the two objects
sobj_cellchat_senmayoCtrl <- subset(Seurat.object,barcodes %in% barcodes_ctrSenmayo)
sobj_cellchat_senmayoSen <- subset(Seurat.object,barcodes %in% barcodes_senSenmayo)

saveRDS(sobj_cellchat_senmayoCtrl,file = "../../out/object/revision/134_sobj_cellchat_senmayoCtrl_VAS.rds")
saveRDS(sobj_cellchat_senmayoSen,file = "../../out/object/revision/134_sobj_cellchat_senmayoSen_VAS.rds")

# SIT
# pull all the cells from CA sample. for the IMM subset pull either the non senescent or the senescent cells
barcodes_senSIT <- LUT_SIT %>%
  filter(pathology_class == "WM_CA") %>%
  mutate(keep = case_when(
    expertAnno.l1 != "VAS" ~ 1,
    expertAnno.l1 == "VAS" & SENEQUANTILE == "YES" ~ 1,
    T~0)) %>%
  filter(keep ==1) %>%
  pull(rowname)

barcodes_ctrSIT <- LUT_SIT %>%
  filter(pathology_class == "WM_CA") %>%
  mutate(keep = case_when(
    expertAnno.l1 != "VAS" ~ 1,
    expertAnno.l1 == "VAS" & SENEQUANTILE == "NO" ~ 1,
    T ~ 0)) %>%
  filter(keep == 1) %>%
  pull(rowname)

length(barcodes_ctrSIT)
length(barcodes_senSIT)

# 3019+182+1010+2864+92+217+6541+570+
# 597
# 194

# save the two objects
sobj_cellchat_SITCtrl <- subset(Seurat.object,barcodes %in% barcodes_ctrSIT)
sobj_cellchat_SITSen <- subset(Seurat.object,barcodes %in% barcodes_senSIT)

saveRDS(sobj_cellchat_SITCtrl,file = "../../out/object/revision/134_sobj_cellchat_SITCtrl_VAS.rds")
saveRDS(sobj_cellchat_SITSen,file = "../../out/object/revision/134_sobj_cellchat_SITSen_VAS.rds")
