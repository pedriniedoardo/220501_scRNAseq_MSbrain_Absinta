# AIM ---------------------------------------------------------------------
# this is the initial step for the counting of senescent cells. In this step I generate the signatures score.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(lemon)

# read in the data --------------------------------------------------------
folder <-  "../../out/table/revision/modules_SENESCENCE_shirmer/"
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

# wrangling ---------------------------------------------------------------
# use wm ctrl as reference
id_signature <- unique(df_modules$signature)
# x <- "senmayo"
list_enumeration <- lapply(id_signature,function(x){
  # try to reference to top10% based on the control dataset
  df_90 <- df_modules %>% 
    filter(signature %in% x,
           pathology %in% "control") %>% 
    group_by(signature,cellid) %>% 
    summarise(thr = quantile(signature_score1,prob=0.90))
  
  # enumerate the cells 090 from the control only samples
  df_en_090_1 <- df_modules %>% 
    select(barcode,signature,pathology,cellid,signature_score1) %>% 
    filter(signature %in% x) %>% 
    left_join(df_90,by = c("cellid","signature")) %>% 
    mutate(pass_thr = signature_score1>thr) %>% 
    group_by(signature,cellid,pathology) %>% 
    summarise(tot_cell = n(),
              tot_sen = sum(pass_thr),.groups = "drop") %>%
    # complete(signature,cellid,pathology, fill = list(tot_sen = 0)) %>%
    # ungroup() %>% 
    mutate(prop_sen = tot_sen/tot_cell)
  
  # count the pro of the contols
  df_en_090_2 <- df_en_090_1 %>% 
    filter(pathology %in% "control") %>% 
    select(signature,cellid,prop_sen_ctrl=prop_sen)
  
  df_en_090_final <- df_en_090_1 %>% 
    left_join(df_en_090_2,by = c("cellid","signature")) %>% 
    mutate(FC = prop_sen/prop_sen_ctrl,
           log2FC = log2(FC))
  
  return(df_en_090_final)
}) %>% 
  setNames(id_signature)

df_en_090_final <- bind_rows(list_enumeration) %>%
  mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active")))
# save the table
df_en_090_final %>% 
  write_tsv("../../out/table/revision/modules_SENESCENCE_shirmer/125_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_senescence_shirmer.tsv")

# plot all the the enumeration
df_en_090_final %>% 
  filter(tot_sen>5) %>% 
  # mutate(origin2 = paste0(disease,"_",origin)) %>% 
  ggplot(aes(y=pathology,x=log2FC))+
  geom_col(aes(fill=tot_sen))+
  # geom_col(aes(width = tot_sen))+
  facet_grid(signature~cellid)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.y.right = element_text(angle = 0))+
  scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
# ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/HVG4000_expertAnno/enumeration_log2Prop_ALL_090_threshold_above5_MSStatus_CXref_cellID_autophagy.pdf"),width = 25,height =15)

# plot a smaller version of the image
# df_en_090_final %>%
#   filter(signature %in% c("CLASSICAL_SASP","Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
#   filter(tot_sen>5) %>%
#   mutate(origin2 = paste0(disease,"_",origin)) %>%
#   ggplot(aes(y=origin2,x=log2FC))+
#   geom_col(aes(fill=tot_sen))+
#   # geom_col(aes(width = tot_sen))+
#   facet_grid(signature~expertAnno.l1)+
#   theme_bw()+
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA),
#         strip.text.y.right = element_text(angle = 0))+
#   scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
# ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/HVG4000_expertAnno/enumeration_log2Prop_ALL_090_threshold_above5_MSStatus_CXref_cellID_small.pdf"),width = 18,height =10)

# reproduce the plot but rename and drop some signature
df_en_090_final %>%
  filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  filter(tot_sen>5) %>%
  # mutate(origin2 = paste0(disease,"_",origin)) %>%
  ggplot(aes(y=pathology,x=log2FC))+
  geom_col(aes(fill=tot_sen))+
  # geom_col(aes(width = tot_sen))+
  facet_grid(signature2~cellid)+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text.y.right = element_text(angle = 0))+
  scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/125_enumeration_log2Prop_ALL_090_threshold_above5_MSStatus_CXref_cellID_small2_shirmer.pdf"),width = 10,height =6)

# -------------------------------------------------------------------------
# run the proportion by cell by sample
# martina recommended to use the cx as reference
folder <-  "../../out/table/revision/modules_SENESCENCE_shirmer/"
# focus on the autophagy only sigantures
pattern_sig <- paste0(names(readRDS("../../data/signatures/senescence_pathways.rds")),collapse = "|")

file <- dir(folder) %>% 
  str_subset(pattern = "Module") %>%
  str_subset(pattern = pattern_sig)

df_modules_cellID <- lapply(file, function(x){
  read_tsv(paste0(folder,x))
}) %>% 
  bind_rows() %>% 
  mutate(disease = factor(disease,levels = c("CTRL","MS"))) %>%
  mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active")))

# LUT_sample <- df_modules_cellID %>% 
#   group_by(sample,pathology,disease) %>% 
#   summarise()

# use wm ctrl as reference
id_signature <- unique(df_modules_cellID$signature)

# x <- "senmayo"
list_enumeration_cellID <- lapply(id_signature,function(x){
  df_90_cellID2 <- df_modules_cellID %>% 
    filter(signature %in% x,
           disease %in% "CTRL") %>% 
    group_by(signature,cellid) %>% 
    summarise(thr = quantile(signature_score1,prob=0.90))
  
  test_summary_cellID <- df_modules_cellID %>% 
    filter(signature %in% x) %>% 
    left_join(df_90_cellID2,c("cellid","signature")) %>% 
    mutate(sen = case_when(signature_score1>thr~1,
                           T~0)) %>% 
    group_by(signature,sample,disease,pathology,cellid,sen) %>% 
    summarise(n = n(),.groups = "drop") %>% 
    # complete(signature,sample,cellid,sen, fill = list(n = 0)) %>%
    ungroup() %>% 
    group_by(signature,sample,disease,pathology,cellid) %>% 
    mutate(tot = sum(n)) %>% 
    ungroup() %>% 
    mutate(prop = n/tot)
  
  # add the sample id
  df_plot <- test_summary_cellID %>% 
    # left_join(LUT_sample,"sample") %>% 
    mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active")))
  
  return(df_plot)
  
}) %>% 
  setNames(id_signature)

enumeration_cellID <- bind_rows(list_enumeration_cellID)

# save the table
enumeration_cellID %>% 
  write_tsv("../../out/table/revision/modules_SENESCENCE_shirmer/125_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_senescence_shirmer.tsv")

# test_summary_cellID %>% 
#   # filter(expertAnno.l1 %in% c(3,5,11)) %>%
#   filter(sen == 1) %>% 
#   ggplot(aes(x=origin,y=prop,col=origin,label=orig.ident)) +
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1)) +
#   geom_text_repel()+
#   facet_grid(expertAnno.l1~disease)+theme_bw()+theme(strip.background = element_blank())+geom_hline(yintercept = 0.1,linetype="dotted",col="black")

# enumeration_cellID %>%
#   # filter(expertAnno.l1 %in% c(3,5,11)) %>%
#   filter(sen == 1) %>% 
#   # ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
#   ggplot(aes(x=pathology_class,y=prop,col=origin)) +
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1)) +
#   # geom_text_repel()+
#   facet_grid(signature~expertAnno.l1,scales = "free")+
#   theme_bw()+
#   theme(strip.background = element_blank(),
#         strip.text.y.right = element_text(angle = 0))+
#   geom_hline(yintercept = 0.1,linetype="dotted",col="black")+theme(axis.text.x = element_text(hjust = 1,angle = 90))
# 
# ggsave(paste0("../../out/image/revision/modules_SENESCENCE/121_enumeration_log2Prop_ALL_090_threshold_above5_MSStatus_CXref_cellID_sampleWise_senescence.pdf"),width = 25,height =15)

# enumeration_cellID %>%
#   # filter(expertAnno.l1 %in% c(3,5,11)) %>%
#   filter(sen == 1) %>% 
#   filter(signature %in% c("CLASSICAL_SASP","Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>% 
#   # ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
#   ggplot(aes(x=pathology_class,y=prop,col=origin)) +
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1)) +
#   # geom_text_repel()+
#   facet_grid(signature~expertAnno.l1,scales = "free")+
#   theme_bw()+
#   theme(strip.background = element_blank(),
#         strip.text.y.right = element_text(angle = 0))+
#   geom_hline(yintercept = 0.1,linetype="dotted",col="black")+theme(axis.text.x = element_text(hjust = 1,angle = 90))
# 
# ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/HVG4000_expertAnno/enumeration_log2Prop_ALL_090_threshold_above5_MSStatus_CXref_cellID_sampleWise_small.pdf"),width = 18,height =10)

# reproduce the plot but rename and drop some signature
enumeration_cellID %>%
  # filter(expertAnno.l1 %in% c(3,5,11)) %>%
  filter(sen == 1) %>%
  filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  # ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
  ggplot(aes(x=pathology,y=prop)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1) +
  # geom_text_repel()+
  facet_grid(signature2~cellid,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 0))+
  geom_hline(yintercept = 0.1,linetype="dotted",col="black")+theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/image/revision/modules_SENESCENCE_shirmer/125_enumeration_log2Prop_ALL_090_threshold_above5_MSStatus_CXref_cellID_sampleWise_small3_shirmer.pdf",width = 8,height =5)

# test account for the zeros ----------------------------------------------
# account for the zeros
# folder <-  "../../out/table/revision/modules_SENESCENCE_shirmer/"
# focus on the autophagy only sigantures
# pattern_sig <- paste0(names(readRDS("../../data/signatures/senescence_pathways.rds")),collapse = "|")

# file <- dir(folder) %>% 
#   str_subset(pattern = "Module") %>%
#   str_subset(pattern = pattern_sig)

# df_modules_cellID <- lapply(file, function(x){
#   read_tsv(paste0(folder,x))
# }) %>% 
#   bind_rows() %>% 
#   mutate(disease = factor(disease,levels = c("CTRL","MS"))) %>%
#   mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active")))

# use wm ctrl as reference
# id_signature <- unique(df_modules_cellID$signature)

LUT_sample <- df_modules_cellID %>%
  group_by(sample,disease,pathology) %>%
  summarise(.groups = "drop")

# x <- "senmayo"
list_enumeration_cellID_zeros <- lapply(id_signature,function(x){
  df_90_cellID2 <- df_modules_cellID %>% 
    filter(signature %in% x,
           disease %in% "CTRL") %>% 
    group_by(signature,cellid) %>% 
    summarise(thr = quantile(signature_score1,prob=0.90))
  
  # tot cells
  tot_cells <- df_modules_cellID %>% 
    filter(signature %in% x) %>% 
    group_by(signature,sample,disease,pathology,cellid) %>% 
    summarise(n = n(),.groups = "drop") %>%
    complete(signature,sample,cellid, fill = list(n = 0)) %>%
    select(-c(disease,pathology)) %>%
    left_join(LUT_sample,by = "sample")
  
  tot_cells_senescence <- df_modules_cellID %>% 
    filter(signature %in% x) %>% 
    left_join(df_90_cellID2,c("cellid","signature")) %>% 
    mutate(sen = case_when(signature_score1>thr~1,
                           T~0)) %>% 
    filter(sen == 1) %>%
    group_by(signature,sample,disease,pathology,cellid) %>% 
    summarise(n_sen = n(),.groups = "drop") %>% 
    complete(signature,sample,cellid, fill = list(n_sen = 0)) %>%
    select(-c(disease,pathology)) %>%
    left_join(LUT_sample,by = "sample")
  
  # join the two tables
  test_summary_cellID <- tot_cells %>%
    left_join(tot_cells_senescence,by = c("signature","sample","cellid","disease","pathology")) %>%
    mutate(prop_sen = n_sen/n)
  
  # add the sample id
  df_plot <- test_summary_cellID %>% 
    # left_join(LUT_sample,"sample") %>% 
    mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active")))
  
  return(df_plot)
  
}) %>% 
  setNames(id_signature)

enumeration_cellID_zeroes <- bind_rows(list_enumeration_cellID_zeros)

# save the table
enumeration_cellID_zeroes %>% 
  write_tsv("../../out/table/revision/modules_SENESCENCE_shirmer/125_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_senescence_shirmer_zeroes.tsv")

# reproduce the plot but rename and drop some signature
enumeration_cellID_zeroes %>%
  # filter(expertAnno.l1 %in% c(3,5,11)) %>%
  filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  # ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
  ggplot(aes(x=pathology,y=prop_sen)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1) +
  # geom_text_repel()+
  facet_grid(signature2~cellid,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 0))+
  geom_hline(yintercept = 0.1,linetype="dotted",col="black")+theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/image/revision/modules_SENESCENCE_shirmer/125_enumeration_log2Prop_ALL_090_threshold_above5_MSStatus_CXref_cellID_sampleWise_small3_shirmer_zeroes.pdf",width = 8,height =5)
