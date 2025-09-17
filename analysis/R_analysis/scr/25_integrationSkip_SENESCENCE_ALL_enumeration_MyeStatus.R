# libraries ---------------------------------------------------------------
library(tidyverse)

# read in the data --------------------------------------------------------
folder <-  "../../out/table/ManualClean/modules_SENESCENCE/"
file <- dir(folder) %>% 
  str_subset(pattern = "Module") %>% 
  str_subset(pattern = "CX_WM_")

df_modules <- lapply(file, function(x){
  read_tsv(paste0(folder,x))
}) %>% 
  bind_rows() %>% 
  mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","WM_Ctrl","WM_NAWM","CX_Mye","WM_CI","WM_CA","CX_Demye","WM_Core")))

# wrangling ---------------------------------------------------------------
# df_en_med_final
# df_modules %>% 
#   filter(signature %in% "senmayo") %>% 
#   group_by(seurat_clusters,pathology) %>% 
#   summarise(n = n())
# df_modules %>% 
#   select(barcode,signature,pathology,seurat_clusters,signature_score1) %>% 
#   filter(signature %in% "senmayo") %>% 
#   left_join(df_med,by = "seurat_clusters") %>% 
#   mutate(pass_thr = signature_score1>thr) %>% 
#   filter(pass_thr==T) %>% 
#   group_by(seurat_clusters,pathology) %>% 
#   summarise(n = n())

# try to reference to top10% based on the control dataset
df_90 <- df_modules %>% 
  filter(signature %in% "senmayo",
         pathology_class %in% "CX_Ctrl") %>% 
  group_by(seurat_clusters) %>% 
  summarise(thr = quantile(signature_score1,prob=0.90))

# enumerate the cells 090 from the control only samples
df_en_090_1 <- df_modules %>% 
  select(barcode,signature,pathology_class,seurat_clusters,signature_score1) %>% 
  filter(signature %in% "senmayo") %>% 
  left_join(df_90,by = "seurat_clusters") %>% 
  mutate(pass_thr = signature_score1>thr) %>% 
  group_by(seurat_clusters,pathology_class) %>% 
  summarise(tot_cell = n(),
            tot_sen = sum(pass_thr)) %>%
  ungroup() %>% 
  mutate(prop_sen = tot_sen/tot_cell)

# count the pro of the contols
df_en_090_2 <- df_en_090_1 %>% 
  filter(pathology_class %in% "CX_Ctrl") %>% 
  select(seurat_clusters,prop_sen_ctrl=prop_sen)

df_en_090_final <- df_en_090_1 %>% 
  left_join(df_en_090_2,by = "seurat_clusters") %>% 
  mutate(FC = prop_sen/prop_sen_ctrl,
         log2FC = log2(FC))

df_en_090_final %>%
  write_tsv("../../out/table/ManualClean/modules_SENESCENCE/enumeration_senmayo_090_CX_WM_PathologyClass_CX_ref.tsv")

# -------------------------------------------------------------------------
# plot using the 0.90 quantile
df_modules %>% 
  filter(signature %in% "senmayo") %>% 
  # mutate(BraakStage=as.factor(BraakStage)) %>% 
  ggplot(aes(x=signature_score1,y=pathology_class,fill=pathology_class))+
  ggridges::geom_density_ridges(alpha=0.5)+
  facet_wrap(~seurat_clusters,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  # scale_fill_manual(values = c("green","yellow","red"))+
  geom_vline(data = df_90,aes(xintercept = thr),linetype="dashed",col="red")
ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/dist_score_ridges_senmayo_090_CX_WM_PathologyClass_CX_ref.pdf"),width = 12,height =9)

# plot the enumeration
df_en_090_final %>% 
  filter(tot_sen>5) %>% 
  ggplot(aes(y=pathology_class,x=log2FC))+
  geom_col(aes(fill=tot_sen))+
  # geom_col(aes(width = tot_sen))+
  facet_wrap(~seurat_clusters)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/enumeration_log2Prop_enumeration_senmayo_090_CX_WM_PathologyClass_CX_ref.pdf"),width = 12,height =9)

# do the same but use WM as reference -------------------------------------
# try to reference to top10% based on the control dataset
df_90_WM <- df_modules %>% 
  filter(signature %in% "senmayo",
         pathology_class %in% "WM_Ctrl") %>% 
  group_by(seurat_clusters) %>% 
  summarise(thr = quantile(signature_score1,prob=0.90))

# enumerate the cells 090 from the control only samples
df_en_090_1_WM <- df_modules %>% 
  select(barcode,signature,pathology_class,seurat_clusters,signature_score1) %>% 
  filter(signature %in% "senmayo") %>% 
  left_join(df_90_WM,by = "seurat_clusters") %>% 
  mutate(pass_thr = signature_score1>thr) %>% 
  group_by(seurat_clusters,pathology_class) %>% 
  summarise(tot_cell = n(),
            tot_sen = sum(pass_thr)) %>%
  ungroup() %>% 
  mutate(prop_sen = tot_sen/tot_cell)

# count the pro of the contols
df_en_090_2_WM <- df_en_090_1_WM %>% 
  filter(pathology_class %in% "WM_Ctrl") %>% 
  select(seurat_clusters,prop_sen_ctrl=prop_sen)

df_en_090_final_WM <- df_en_090_1_WM %>% 
  left_join(df_en_090_2_WM,by = "seurat_clusters") %>% 
  mutate(FC = prop_sen/prop_sen_ctrl,
         log2FC = log2(FC))

df_en_090_final_WM %>%
  write_tsv("../../out/table/ManualClean/modules_SENESCENCE/enumeration_senmayo_090_CX_WM_PathologyClass_WM_ref.tsv")

# -------------------------------------------------------------------------
# plot using the 0.90 quantile
df_modules %>% 
  filter(signature %in% "senmayo") %>% 
  # mutate(BraakStage=as.factor(BraakStage)) %>% 
  ggplot(aes(x=signature_score1,y=pathology_class,fill=pathology_class))+
  ggridges::geom_density_ridges(alpha=0.5)+
  facet_wrap(~seurat_clusters,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  # scale_fill_manual(values = c("green","yellow","red"))+
  geom_vline(data = df_90_WM,aes(xintercept = thr),linetype="dashed",col="red")
ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/dist_score_ridges_senmayo_090_CX_WM_PathologyClass_WM_ref.pdf"),width = 12,height =9)

# plot the enumeration
df_en_090_final_WM %>% 
  filter(tot_sen>5) %>% 
  ggplot(aes(y=pathology_class,x=log2FC))+
  geom_col(aes(fill=tot_sen))+
  # geom_col(aes(width = tot_sen))+
  facet_wrap(~seurat_clusters)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/enumeration_log2Prop_enumeration_senmayo_090_CX_WM_PathologyClass_WM_ref.pdf"),width = 12,height =9)
