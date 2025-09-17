# libraries ---------------------------------------------------------------
library(tidyverse)
library(ggrepel)

# read in the data --------------------------------------------------------
folder <-  "../../out/table/ManualClean/modules_SENESCENCE/HVG4000/"
file <- dir(folder) %>% 
  str_subset(pattern = "Module")

df_modules <- lapply(file, function(x){
  read_tsv(paste0(folder,x))
}) %>% 
  bind_rows() %>% 
  mutate(disease = factor(disease,levels = c("CTRL","MS")))

# wrangling ---------------------------------------------------------------
# use wm ctrl as reference

# try to reference to top10% based on the control dataset
df_90 <- df_modules %>% 
  filter(signature %in% "senmayo",
         disease %in% "CTRL",
         origin == "wm") %>% 
  group_by(seurat_clusters) %>% 
  summarise(thr = quantile(signature_score1,prob=0.90))

# enumerate the cells 090 from the control only samples
df_en_090_1 <- df_modules %>% 
  select(barcode,signature,disease,seurat_clusters,signature_score1,origin) %>% 
  filter(signature %in% "senmayo") %>% 
  left_join(df_90,by = "seurat_clusters") %>% 
  mutate(pass_thr = signature_score1>thr) %>% 
  group_by(seurat_clusters,disease,origin) %>% 
  summarise(tot_cell = n(),
            tot_sen = sum(pass_thr)) %>%
  ungroup() %>% 
  mutate(prop_sen = tot_sen/tot_cell)

# count the pro of the contols
df_en_090_2 <- df_en_090_1 %>% 
  filter(disease %in% "CTRL",
         origin == "wm") %>% 
  select(seurat_clusters,prop_sen_ctrl=prop_sen)

df_en_090_final <- df_en_090_1 %>% 
  left_join(df_en_090_2,by = "seurat_clusters") %>% 
  mutate(FC = prop_sen/prop_sen_ctrl,
         log2FC = log2(FC))

df_en_090_final %>%
  write_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refWM.tsv")

# read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refWM.tsv")

# -------------------------------------------------------------------------
# plot using the 0.90 quantile
df_modules %>% 
  filter(signature %in% "senmayo") %>% 
  # mutate(BraakStage=as.factor(BraakStage)) %>% 
  ggplot(aes(x=signature_score1,y=paste0(disease,"_",origin),fill=disease))+
  ggridges::geom_density_ridges(alpha=0.5)+
  facet_wrap(~seurat_clusters,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_fill_manual(values = c("green","red"))+
  geom_vline(data = df_90,aes(xintercept = thr),linetype="dashed",col="red")
ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/dist_score_ridges_senmayo_090_threshold_MSStatus.pdf"),width = 12,height =9)

# plot the enumeration
df_en_090_final %>% 
  filter(tot_sen>5) %>% 
  ggplot(aes(y=paste0(disease,"_",origin),x=log2FC))+
  geom_col(aes(fill=tot_sen))+
  # geom_col(aes(width = tot_sen))+
  facet_wrap(~seurat_clusters)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/enumeration_log2Prop_senmayo_090_threshold_above5_MSStatus.pdf"),width = 12,height =9)

# -------------------------------------------------------------------------
# use cx ctrl as reference

# try to reference to top10% based on the control dataset
df_90 <- df_modules %>% 
  filter(signature %in% "senmayo",
         disease %in% "CTRL",
         origin == "cortex") %>% 
  group_by(seurat_clusters) %>% 
  summarise(thr = quantile(signature_score1,prob=0.90))

# enumerate the cells 090 from the control only samples
df_en_090_1 <- df_modules %>% 
  select(barcode,signature,disease,seurat_clusters,signature_score1,origin) %>% 
  filter(signature %in% "senmayo") %>% 
  left_join(df_90,by = "seurat_clusters") %>% 
  mutate(pass_thr = signature_score1>thr) %>% 
  group_by(seurat_clusters,disease,origin) %>% 
  summarise(tot_cell = n(),
            tot_sen = sum(pass_thr)) %>%
  ungroup() %>% 
  mutate(prop_sen = tot_sen/tot_cell)

# count the pro of the contols
df_en_090_2 <- df_en_090_1 %>% 
  filter(disease %in% "CTRL",
         origin == "cortex") %>% 
  select(seurat_clusters,prop_sen_ctrl=prop_sen)

df_en_090_final <- df_en_090_1 %>% 
  left_join(df_en_090_2,by = "seurat_clusters") %>% 
  mutate(FC = prop_sen/prop_sen_ctrl,
         log2FC = log2(FC))

df_en_090_final %>%
  write_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refCX.tsv")

read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refCX.tsv")

# -------------------------------------------------------------------------
# plot using the 0.90 quantile
df_modules %>% 
  filter(signature %in% "senmayo") %>% 
  # mutate(BraakStage=as.factor(BraakStage)) %>% 
  ggplot(aes(x=signature_score1,y=paste0(disease,"_",origin),fill=disease))+
  ggridges::geom_density_ridges(alpha=0.5)+
  facet_wrap(~seurat_clusters,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_fill_manual(values = c("green","red"))+
  geom_vline(data = df_90,aes(xintercept = thr),linetype="dashed",col="red")
ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/dist_score_ridges_senmayo_090_threshold_MSStatus.pdf"),width = 12,height =9)

# plot the enumeration
df_en_090_final %>% 
  filter(tot_sen>5) %>% 
  ggplot(aes(y=paste0(disease,"_",origin),x=log2FC))+
  geom_col(aes(fill=tot_sen))+
  # geom_col(aes(width = tot_sen))+
  facet_wrap(~seurat_clusters)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/enumeration_log2Prop_senmayo_090_threshold_above5_MSStatus.pdf"),width = 12,height =9)

# load the ratio and quantifications --------------------------------------
# compare total proportions
read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refWM.tsv") %>% 
  filter(disease == "MS") %>% 
  
  ggplot(aes(x=factor(seurat_clusters),y=prop_sen,fill=origin))+geom_col(position = "dodge")+facet_grid(~disease)+theme_bw()+theme(strip.background = element_blank())+geom_hline(yintercept = 0.1,linetype="dotted",col="black")

read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refCX.tsv") %>% 
  filter(disease == "MS") %>% 
  
  ggplot(aes(x=factor(seurat_clusters),y=prop_sen,fill=origin))+geom_col(position = "dodge")+facet_grid(~disease)+theme_bw()+theme(strip.background = element_blank())+geom_hline(yintercept = 0.1,linetype="dotted",col="black")


# pull the senescent cells from each cluster ------------------------------
# pull the senescent cells, how many are there from individula samples
df_90 <- df_modules %>% 
  filter(signature %in% "senmayo",
         disease %in% "CTRL",
         origin == "wm") %>% 
  group_by(seurat_clusters) %>% 
  summarise(thr = quantile(signature_score1,prob=0.90))

test_summary <- df_modules %>%
  filter(signature %in% "senmayo") %>% 
  left_join(df_90,"seurat_clusters") %>% 
  mutate(sen = case_when(signature_score1>thr~1,
         T~0)) %>% 
  group_by(orig.ident,origin,disease,seurat_clusters,sen) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(orig.ident,origin,disease,seurat_clusters) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop = n/tot)

test_summary %>% 
  filter(seurat_clusters %in% c(3,5,11)) %>%
  filter(sen == 1) %>% 
  ggplot(aes(x=origin,y=prop,col=origin)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1)) +
  facet_grid(seurat_clusters~disease)+theme_bw()+theme(strip.background = element_blank())+geom_hline(yintercept = 0.1,linetype="dotted",col="black")

test_summary %>% 
  filter(seurat_clusters %in% c(3,5,11)) %>%
  filter(sen == 1) %>% 
  ggplot(aes(x=origin,y=prop,col=origin,label=orig.ident)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1)) +
  geom_text_repel()+
  facet_grid(seurat_clusters~disease)+theme_bw()+theme(strip.background = element_blank())+geom_hline(yintercept = 0.1,linetype="dotted",col="black")

# -------------------------------------------------------------------------
df_90 <- df_modules %>% 
  filter(signature %in% "senmayo",
         disease %in% "CTRL",
         origin == "cortex") %>% 
  group_by(seurat_clusters) %>% 
  summarise(thr = quantile(signature_score1,prob=0.90))

test_summary2 <- df_modules %>% 
  filter(signature %in% "senmayo") %>% 
  left_join(df_90,"seurat_clusters") %>% 
  mutate(sen = case_when(signature_score1>thr~1,
                         T~0)) %>% 
  group_by(orig.ident,origin,disease,seurat_clusters,sen) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(orig.ident,origin,disease,seurat_clusters) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop = n/tot)

test_summary2 %>% 
  filter(seurat_clusters %in% c(3,5,11)) %>%
  filter(sen == 1) %>% 
  ggplot(aes(x=origin,y=prop,col=origin,label=orig.ident)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1)) +
  geom_text_repel()+
  facet_grid(seurat_clusters~disease)+theme_bw()+theme(strip.background = element_blank())+geom_hline(yintercept = 0.1,linetype="dotted",col="black")

# tailored plots ----------------------------------------------------------
# martina recommended to use the cx as reference
folder <-  "../../out/table/ManualClean/modules_SENESCENCE/HVG4000/"
file <- dir(folder) %>% 
  str_subset(pattern = "Module")

df_modules_cellID <- lapply(file, function(x){
  read_tsv(paste0(folder,x))
}) %>% 
  bind_rows() %>% 
  mutate(disease = factor(disease,levels = c("CTRL","MS"))) %>% 
  mutate(cell_id = case_when(seurat_clusters %in% c(0,4,14)~"OLIGO",
                             seurat_clusters %in% c(9)~"OPC",
                             seurat_clusters %in% c(3,12)~"ASTRO",
                             seurat_clusters %in% c(5)~"IMMUNE",
                             seurat_clusters %in% c(13)~"LYM",
                             seurat_clusters %in% c(11)~"VAS",
                             seurat_clusters %in% c(1, 2, 10,6)~"EXC NEU",
                             seurat_clusters %in% c(7,8)~"INH NEU"))


# try to reference to top10% based on the control dataset
df_90_cellID <- df_modules_cellID %>% 
  filter(signature %in% "senmayo",
         disease %in% "CTRL",
         origin == "cortex") %>% 
  group_by(cell_id) %>% 
  summarise(thr = quantile(signature_score1,prob=0.90))

# enumerate the cells 090 from the control only samples
df_en_090_1_cellID <- df_modules_cellID %>% 
  select(barcode,signature,disease,cell_id,signature_score1,origin) %>% 
  filter(signature %in% "senmayo") %>% 
  left_join(df_90_cellID,by = "cell_id") %>% 
  mutate(pass_thr = signature_score1>thr) %>% 
  group_by(cell_id,disease,origin,signature) %>% 
  summarise(tot_cell = n(),
            tot_sen = sum(pass_thr)) %>%
  ungroup() %>% 
  mutate(prop_sen = tot_sen/tot_cell)

# count the pro of the contols
df_en_090_2_cellID <- df_en_090_1_cellID %>% 
  filter(disease %in% "CTRL",
         origin == "cortex") %>% 
  select(cell_id,prop_sen_ctrl=prop_sen)

df_en_090_final_cellID <- df_en_090_1_cellID %>% 
  left_join(df_en_090_2_cellID,by = "cell_id") %>% 
  mutate(FC = prop_sen/prop_sen_ctrl,
         log2FC = log2(FC))

df_en_090_final_cellID %>%
  write_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refCX_cellID.tsv")

# -------------------------------------------------------------------------
# plot using the 0.90 quantile
df_modules_cellID %>% 
  filter(signature %in% "senmayo") %>% 
  # mutate(BraakStage=as.factor(BraakStage)) %>% 
  ggplot(aes(x=signature_score1,y=paste0(disease,"_",origin),fill=disease))+
  ggridges::geom_density_ridges(alpha=0.5)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_fill_manual(values = c("green","red"))+
  geom_vline(data = df_90_cellID,aes(xintercept = thr),linetype="dashed",col="red")
ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/dist_score_ridges_senmayo_090_threshold_MSStatus_cellID.pdf"),width = 12,height =9)

# plot the enumeration
df_en_090_final_cellID %>% 
  filter(tot_sen>5) %>% 
  ggplot(aes(y=paste0(disease,"_",origin),x=log2FC))+
  geom_col(aes(fill=tot_sen))+
  # geom_col(aes(width = tot_sen))+
  facet_wrap(~cell_id)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/enumeration_log2Prop_senmayo_090_threshold_above5_MSStatus_cellID.pdf"),width = 12,height =9)

# load the ratio and quantifications --------------------------------------
# compare total proportions
read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refCX_cellID.tsv") %>% 
  filter(disease == "MS") %>% 
  
  ggplot(aes(x=factor(cell_id),y=prop_sen,fill=origin))+geom_col(position = "dodge")+facet_grid(~disease)+theme_bw()+theme(strip.background = element_blank())+geom_hline(yintercept = 0.1,linetype="dotted",col="black")

# pull the senescent cells from each cluster ------------------------------
# pull the senescent cells, how many are there from individula samples
df_90_cellID2 <- df_modules_cellID %>% 
  filter(signature %in% "senmayo",
         disease %in% "CTRL",
         origin == "cortex") %>% 
  group_by(cell_id) %>% 
  summarise(thr = quantile(signature_score1,prob=0.90))

# df_modules_cellID %>% 
#   filter(signature %in% "senmayo",
#          disease %in% "CTRL",
#          origin == "cortex") %>%
#   filter(cell_id == "LYM") %>% 
#   select(signature_score1) %>% 
#   summarise(thr = quantile(signature_score1,prob=0.90))

# df_modules_cellID %>% 
#   filter(signature_score1>0.0659,cell_id=="LYM",origin=="cortex")

test_summary_cellID <- df_modules_cellID %>% 
  filter(signature %in% "senmayo") %>% 
  left_join(df_90_cellID2,"cell_id") %>% 
  mutate(sen = case_when(signature_score1>thr~1,
                         T~0)) %>% 
  group_by(signature,orig.ident,origin,disease,cell_id,sen) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(signature,orig.ident,origin,disease,cell_id) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop = n/tot)

test_summary_cellID %>%
  write_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refCX_cellID_sampleWise.tsv")

# test_summary_cellID %>% 
#   filter(cell_id == "LYM") %>% 
#   group_by(origin,disease,cell_id,sen) %>% 
#   summarise(n = sum(n)) %>% 
#   ungroup() %>% 
#   group_by(origin,disease,cell_id) %>% 
#   mutate(tot = sum(n)) %>% 
#   ungroup() %>% 
#   mutate(prop = n/tot)
# 
# df_modules_cellID %>% 
#   filter(origin=="cortex",disease=="CTRL",cell_id=="LYM")
# 
# test_summary_cellID %>% 
#   filter(origin=="cortex",disease=="CTRL",cell_id=="LYM")
# 
# df_modules_cellID %>% 
#   left_join(df_90_cellID2,"cell_id") %>% 
#   mutate(sen = case_when(signature_score1>thr~1,
#                          T~0)) %>% 
#   group_by(orig.ident,origin,disease,cell_id,sen) %>% 
#   summarise(n = n()) %>% 
#   filter(cell_id == "LYM")
# 
# test_summary_cellID %>% 
#   filter(cell_id == "LYM",origin=="cortex",disease=="CTRL")

test_summary_cellID %>% 
  # filter(seurat_clusters %in% c(3,5,11)) %>%
  filter(sen == 1) %>% 
  ggplot(aes(x=origin,y=prop,col=origin)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1)) +
  facet_grid(cell_id~disease)+theme_bw()+theme(strip.background = element_blank())+geom_hline(yintercept = 0.1,linetype="dotted",col="black")

test_summary_cellID %>% 
  # filter(cell_id %in% c(3,5,11)) %>%
  filter(sen == 1) %>% 
  ggplot(aes(x=origin,y=prop,col=origin,label=orig.ident)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1)) +
  geom_text_repel()+
  facet_grid(cell_id~disease)+theme_bw()+theme(strip.background = element_blank())+geom_hline(yintercept = 0.1,linetype="dotted",col="black")

# add the sample id
LUT_sample <- read_tsv("../../out/table/ManualClean/meta_data.combined_WM_CX_harmonySkipIntegAllSoupX.tsv") %>% 
  group_by(orig.ident,pathology_class) %>% 
  summarise()

test_summary_cellID %>% 
  left_join(LUT_sample,"orig.ident") %>% 
  mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Demye","CX_Mye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core"))) %>% 
  # filter(cell_id %in% c(3,5,11)) %>%
  filter(sen == 1) %>% 
  ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1)) +
  #geom_text_repel()+
  facet_grid(cell_id~origin,scales = "free")+theme_bw()+theme(strip.background = element_blank())+geom_hline(yintercept = 0.1,linetype="dotted",col="black")+theme(axis.text.x = element_text(hjust = 1,angle = 90))

test_summary_cellID %>% 
  left_join(LUT_sample,"orig.ident") %>% 
  mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Demye","CX_Mye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core"))) %>% 
  write_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refCX_cellID_sampleWise2.tsv")
