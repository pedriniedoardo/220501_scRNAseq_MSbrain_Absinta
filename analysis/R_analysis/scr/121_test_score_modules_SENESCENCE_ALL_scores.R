# AIM ---------------------------------------------------------------------
# Martina asked to try to provide the average signature score per cell type per sample

# libraries ---------------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(lemon)
library(Seurat)

# read in the data --------------------------------------------------------
folder <-  "../../out/table/revision/modules_SENESCENCE/"
# focus on the autophagy only sigantures
pattern_sig <- paste0(names(readRDS("../../data/signatures/senescence_pathways.rds")),collapse = "|")

file <- dir(folder) %>% 
  str_subset(pattern = "Module") %>%
  str_subset(pattern = pattern_sig)

# they already have the expert annotation suggested by martina
df_modules <- lapply(file, function(x){
  read_tsv(paste0(folder,x))
}) %>% 
  bind_rows() %>% 
  mutate(disease = factor(disease,levels = c("CTRL","MS")))

# read in the dataset
data.combined <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")

# save the plot
DimPlot(data.combined, group.by = "expertAnno.l1",raster=T,label = T)
dim(data.combined)

# merge the data
# add also the coordinates for the UMAP
meta_full <- data.combined@meta.data %>%
  rownames_to_column("rowname")

# summarise the full metadata
dim(meta_full)
dim(data.combined)

# load the sample metadata
LUT_sample <- read_tsv("../../out/table/revision/120_meta.data_AnnotationSCType_manualAnnotation.tsv") %>% 
  group_by(orig.ident,pathology_class) %>% 
  summarise()

# wrangling ---------------------------------------------------------------
# calculate the average score per signature per cell type per sample
df_avg_score <- df_modules %>%
  group_by(signature,origin,pathology_class,orig.ident,disease,expertAnno.l1) %>%
  summarise(avg_score = mean(signature_score1))

# check if average score and proportion are correlated
enumeration_cellID_full <- read_tsv("../../out/table/revision/modules_SENESCENCE/121_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_senescence.tsv") %>%
  filter(sen==1)

left_join(df_avg_score,enumeration_cellID_full,by = c("signature","origin","pathology_class","orig.ident","disease","expertAnno.l1")) %>%
  filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  ggplot(aes(x=prop,y=avg_score))+geom_smooth(method = "lm")+geom_point()+facet_wrap(expertAnno.l1~signature,scales = "free",ncol = 4)+theme_bw()+theme(strip.background = element_blank())
ggsave("../../out/image/121_scatter_SENMAYO_prop_vs_score.pdf",height = 20,width = 10)

# plot all the the scores
df_avg_score %>%
  # filter(expertAnno.l1 %in% c(3,5,11)) %>%
  filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  # ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
  ggplot(aes(x=pathology_class,y=avg_score,col=origin)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1)) +
  # geom_text_repel()+
  facet_grid(signature2~expertAnno.l1,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 0))+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave(paste0("../../out/image/121_score_cellID_sampleWise_small3.pdf"),width = 18,height =8)

# export the table with proportions with all metadata
df_share <- df_avg_score %>%
  filter(signature %in% c("senmayo")) %>%
  left_join(meta_full %>%
              group_by(orig.ident,sample_id,facility,origin,disease,sex,age,pathology,pathology_class,patient,plaque,PMI,origin_alt) %>%
              summarise(),by = c("orig.ident","origin"))

df_share %>%
  write_tsv("../../out/table/revision/121_SENMAYO_score_table.tsv")  
