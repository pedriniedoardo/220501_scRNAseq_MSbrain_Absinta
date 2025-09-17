# AIM ---------------------------------------------------------------------
# this is the initial step for the counting of senescent cells. In this step I generate the signatures score.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(lemon)

# read in the data --------------------------------------------------------
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

# add the young vs old label based on the classification defined by martina
df_modules_cellID_fix <- 
  df_modules_cellID %>%
  # focus only on the control samples
  filter(disease == "CTRL") %>%
  # define young and old
  mutate(age_cat = case_when(age < 40 ~ "young",
                             age > 65 ~ "old")) %>%
  # remove all the non labelled samples
  filter(!is.na(age_cat))

# LUT_sample <- df_modules_cellID %>% 
#   group_by(sample,pathology,disease) %>% 
#   summarise()

# use wm ctrl as reference
id_signature <- unique(df_modules_cellID$signature)

# x <- "senmayo"
list_enumeration_cellID <- lapply(id_signature,function(x){
  df_90_cellID2 <- df_modules_cellID_fix %>% 
    filter(signature %in% x,
           age_cat %in% c("young")) %>% 
    group_by(signature,cellid) %>% 
    summarise(thr = quantile(signature_score1,prob=0.90))
  
  test_summary_cellID <- df_modules_cellID_fix %>% 
    filter(signature %in% x) %>% 
    left_join(df_90_cellID2,c("cellid","signature")) %>% 
    mutate(sen = case_when(signature_score1>thr~1,
                           T~0)) %>% 
    group_by(signature,sample,disease,pathology,age_cat,cellid,sen) %>% 
    summarise(n = n()) %>% 
    ungroup() %>% 
    group_by(signature,sample,disease,pathology,age_cat,cellid) %>% 
    mutate(tot = sum(n)) %>% 
    ungroup() %>% 
    mutate(prop = n/tot)
  
  # add the sample id
  df_plot <- test_summary_cellID %>% 
    # left_join(LUT_sample,"sample") %>% 
    mutate(age_cat = factor(age_cat,levels = c("young","old")))
  
  return(df_plot)
  
}) %>% 
  setNames(id_signature)

enumeration_cellID <- bind_rows(list_enumeration_cellID)

# save the table
enumeration_cellID %>% 
  write_tsv("../../out/table/revision/modules_SENESCENCE_shirmer/125_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_senescence_shirmer_young_vs_old.tsv")


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
  ggplot(aes(x=age_cat,y=prop)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1) +
  # geom_text_repel()+
  facet_grid(signature2~cellid,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 0))+
  geom_hline(yintercept = 0.1,linetype="dotted",col="black")+theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/image/revision/modules_SENESCENCE_shirmer/125_enumeration_log2Prop_ALL_090_threshold_above5_MSStatus_CXref_cellID_sampleWise_small3_shirmer_young_vs_old.pdf",width = 8,height =5)
