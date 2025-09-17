# AIM ---------------------------------------------------------------------
# use a similar approach for the SENESCENECE enumeration, but focus on the control samples alone. in particular focus on the young vs old samples

# libraries ---------------------------------------------------------------
library(tidyverse)
library(broom)

# read in the data --------------------------------------------------------
# folder with the quantification
folder <-  "../../out/table/revision/modules_SENESCENCE_shirmer/"
# folder from absinta data
folder2 <-  "../../out/table/revision/modules_SENESCENCE/"

# focus on the autophagy only sigantures
pattern_sig <- paste0(names(readRDS("../../data/signatures/senescence_pathways.rds")),collapse = "|")

file <- dir(folder) %>% 
  str_subset(pattern = "Module") %>%
  str_subset(pattern = pattern_sig)

file2 <- dir(folder2) %>% 
  str_subset(pattern = "Module") %>%
  str_subset(pattern = pattern_sig)

# shirmer data
df_modules_cellID <- lapply(file, function(x){
  read_tsv(paste0(folder,x))
}) %>% 
  bind_rows() %>% 
  mutate(disease = factor(disease,levels = c("CTRL","MS"))) %>%
  mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active")))

# absinta data
df_modules_cellID2 <- lapply(file2, function(x){
  read_tsv(paste0(folder2,x))
}) %>% 
  bind_rows() %>% 
  mutate(disease = factor(disease,levels = c("CTRL","MS"))) %>%
  mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Demye","CX_Mye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core")))

# wrangling ---------------------------------------------------------------
# approach 1 check the correlation between age and senmayo score per cell type
# in theory we would expect a positive correlation between the score and the age in the control samples
# average the score per sample and per cellid
df_modules_cellID %>%
  filter(disease == "CTRL") %>%
  group_by(sample,age,sex) %>%
  summarise() %>%
  arrange(age)
df_sample_cellid <- df_modules_cellID %>%
  filter(disease == "CTRL") %>%
  group_by(sample,cellid,age,sex,signature) %>%
  summarise(avg_score = mean(signature_score1))

df_modules_cellID2 %>%
  filter(disease == "CTRL") %>%
  group_by(orig.ident,age,sex,origin) %>%
  summarise() %>%
  arrange(age)
df_sample_cellid2 <- df_modules_cellID2 %>%
  filter(disease == "CTRL") %>%
  group_by(orig.ident,origin,expertAnno.l1,age,sex,signature) %>%
  summarise(avg_score = mean(signature_score1))

# plot the data
# shirmer
df_sample_cellid %>%
  filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  ggplot(aes(x=age,y=avg_score))+theme_bw()+
  facet_wrap(cellid~signature2,scales = "free",ncol=4)+geom_smooth(method = "lm")+geom_point(shape=1)+
  theme(strip.background = element_blank())
ggsave("../../out/image/revision/125_correlation_signature_age_shirmer.pdf",width = 10,height = 20)

# absinta
df_sample_cellid2 %>%
  filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  ggplot(aes(x=age,y=avg_score))+theme_bw()+
  facet_wrap(expertAnno.l1~signature2,scales = "free",ncol=4)+geom_smooth(method = "lm",se=T)+geom_point(shape=1)+
  theme(strip.background = element_blank())
ggsave("../../out/image/revision/125_correlation_signature_age_absinta.pdf",width = 10,height = 22.5)

# martina suggeseted to plot only the senmayo column in two rows
df_sample_cellid2 %>%
  filter(signature %in% c("senmayo")) %>%
  filter(!expertAnno.l1 %in% c("EPENDYMA")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  ggplot(aes(x=age,y=avg_score))+theme_bw()+
  facet_wrap(expertAnno.l1~signature2,scales = "free",ncol=4)+geom_smooth(method = "lm",se=T)+geom_point(shape=1)+
  theme(strip.background = element_blank())
ggsave("../../out/image/revision/125_correlation_signature_age_absinta_test.pdf",width = 12,height = 6)

# run the test per cell type
df_sample_cellid2 %>%
  filter(signature %in% c("senmayo")) %>%
  filter(!expertAnno.l1 %in% c("EPENDYMA")) %>%
  group_by(expertAnno.l1) %>%
  nest() %>%
  mutate(corr = map(data,function(x){
    df_out <- cor.test(x$age,x$avg_score)
    return(tidy(df_out))
  })) %>%
  unnest(corr)

# try to plot them all in one scatter
df_sample_cellid2 %>%
  filter(signature %in% c("senmayo")) %>%
  filter(!expertAnno.l1 %in% c("EPENDYMA")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  ggplot(aes(x=age,y=avg_score))+theme_bw()+
  # facet_wrap(expertAnno.l1~signature2,scales = "free",ncol=4)+
  geom_smooth(method = "lm",se=T)+geom_point(shape=1,position = position_jitter(width = 0.5))+
  # geom_smooth(method = "lm",se=T)+geom_point(shape=1)+
  theme(strip.background = element_blank())
ggsave("../../out/image/revision/125_correlation_signature_age_absinta_test_global.pdf",width = 4,height = 4)

df_corr_test <- df_sample_cellid2 %>%
  filter(signature %in% c("senmayo")) %>%
  filter(!expertAnno.l1 %in% c("EPENDYMA"))

cor.test(df_corr_test$age,df_corr_test$avg_score)

test_res <- cor.test(df_corr_test$age,df_corr_test$avg_score)
tidy(test_res)

# absinta with origin
df_sample_cellid2 %>%
  filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  ggplot(aes(x=age,y=avg_score,col=origin))+theme_bw()+
  facet_wrap(expertAnno.l1~signature2,scales = "free",ncol=4)+geom_smooth(method = "lm",aes(col=origin),se=F)+geom_point(shape=1)+
  theme(strip.background = element_blank())
ggsave("../../out/image/revision/125_correlation_signature_age_absinta_origin.pdf",width = 10,height = 22.5)

# -------------------------------------------------------------------------
# try to put together in a global plot, both absinta and shirmer
df_corr_tot <- bind_rows(
  df_sample_cellid2 %>%
    mutate(dataset = "Absinta") %>%
    ungroup() %>%
    select(sample = orig.ident, cellid = expertAnno.l1,dataset,age,avg_score,signature),
  df_sample_cellid %>%
    mutate(dataset = "Shirmer") %>%
    ungroup() %>%
    mutate(sample = as.character(sample)) %>%
    select(sample, cellid,dataset,age,avg_score,signature)
  )

p1 <- df_corr_tot %>%
  filter(signature %in% c("senmayo")) %>%
  filter(!cellid %in% c("EPENDYMA")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  ggplot(aes(x=age,y=avg_score)) + theme_bw() +
  # facet_wrap(expertAnno.l1~signature2,scales = "free",ncol=4)+
  geom_smooth(method = "lm",se=T)+geom_point(shape=1,position = position_jitter(width = 0.5)) +
  # geom_smooth(method = "lm",se=T)+geom_point(shape=1)+
  theme(strip.background = element_blank()) +
  ggtitle("Global correlation") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- df_corr_tot %>%
  filter(signature %in% c("senmayo")) %>%
  filter(!cellid %in% c("EPENDYMA")) %>%
  mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  ggplot(aes(x=age,y=avg_score,col=dataset)) + theme_bw() +
  # facet_wrap(expertAnno.l1~signature2,scales = "free",ncol=4)+
  geom_smooth(method = "lm",se=T)+geom_point(shape=1,position = position_jitter(width = 0.5)) +
  # geom_smooth(method = "lm",se=T)+geom_point(shape=1)+
  theme(strip.background = element_blank())+
  ggtitle("Global correlation by dataset") +
  theme(plot.title = element_text(hjust = 0.5))

p1+p2

df_corr_test <- df_sample_cellid2 %>%
  filter(signature %in% c("senmayo")) %>%
  filter(!expertAnno.l1 %in% c("EPENDYMA"))

cor.test(df_corr_test$age,df_corr_test$avg_score)
