# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)

# read the data -----------------------------------------------------------
# read in the dataset
data.combined <- readRDS("../../data/schirmer.rds")
DimPlot(data.combined,raster = T,reduction = "curated",label = T)
DefaultAssay(data.combined) <- "RNA"
dim(data.combined)

# pull the annotation derived from the SIT algorithm
df_SIT <- read_tsv("../../out/table/revision/124_meta_SIT_WM_CX_harmonySkipIntegration_AllSoupX_test_shirmer.tsv")

# merge the data
meta_full <- data.combined@meta.data %>%
  rownames_to_column("rowname") %>%
  left_join(df_SIT %>%
              dplyr::select("rowname",
                            "KEGG_CELL_CYCLE","Cell_cycle_arrest_signatureK","REACTOME_CELL_CYCLE","Cell_cycle_arrest_signatureR","WP_CELL_CYCLE","Cell_cycle_arrest_signatureWP",
                            "Senescence_signature1",
                            "Senescence_scoreK",
                            "Senescence_scoreWP",
                            "Senescence_scoreR",
                            "Senescence_scoreKQUANTILE",
                            "Senescence_scoreRQUANTILE",
                            "Senescence_scoreWPQUANTILE",
                            "SENEQUANTILE","cellid","disease"),by = "rowname")

# summarise the full metadata
dim(meta_full)
dim(data.combined)

# use the sema criterion used with the senmayo approach
# test_summary_cellID <- df_modules_cellID %>% 
#   filter(signature %in% x) %>% 
#   left_join(df_90_cellID2,c("expertAnno.l1","signature")) %>% 
#   mutate(sen = case_when(signature_score1>thr~1,
#                          T~0)) %>% 
#   group_by(signature,orig.ident,origin,disease,expertAnno.l1,sen) %>% 
#   summarise(n = n()) %>% 
#   ungroup() %>% 
#   group_by(signature,orig.ident,origin,disease,expertAnno.l1) %>% 
#   mutate(tot = sum(n)) %>% 
#   ungroup() %>% 
#   mutate(prop = n/tot)

enumeration_cellID <- meta_full %>%
  group_by(sample,pathology,cellid,SENEQUANTILE) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(sample,pathology,cellid) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot) %>%
  mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active"))) %>%
  mutate(signature ="SIT")

enumeration_cellID %>% 
  write_tsv("../../out/table/revision/124_enumeration_SIT_WM_CX_harmonySkipIntegAllSoupX_expertAnno_shirmer.tsv")

# plot the data to show the proportions of senescent cells per condition
enumeration_cellID %>%
  # filter(expertAnno.l1 %in% c(3,5,11)) %>%
  filter(SENEQUANTILE == "YES") %>% 
  # filter(signature %in% c("CLASSICAL_SASP","Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>% 
  # ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
  ggplot(aes(x=pathology,y=prop)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1)) +
  # geom_text_repel()+
  # facet_grid(signature~expertAnno.l1,scales = "free")+
  facet_wrap(signature~cellid,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 0))+
  # geom_hline(yintercept = 0.1,linetype="dotted",col="black")+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))

ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/124_enumeration_SIT_WM_CX_harmonySkipIntegAllSoupX_expertAnno_shirmer.pdf"),width = 9,height =8)

# check the correlation between SIT and senmayo ---------------------------
# join the two datasets pull only the senescence cells from both
df_senmayo <- read_tsv("../../out/table/revision/modules_SENESCENCE_shirmer/125_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_senescence_shirmer.tsv") %>%
  filter(signature == "senmayo") %>%
  filter(sen == 1)


df_SIT <- read_tsv("../../out/table/revision/124_enumeration_SIT_WM_CX_harmonySkipIntegAllSoupX_expertAnno_shirmer.tsv") %>%
  filter(SENEQUANTILE == "YES")

# check the reason for the discrepancies in the two dataset. the call for senescence is not the same in both, therefore it might be absent in one estimate compard to the other
full_join(df_senmayo %>%
            group_by(sample,pathology) %>%
            summarise(n = n()),
          df_SIT %>%
            group_by(sample,pathology) %>%
            summarise(n = n()),by = c("sample","pathology"),suffix = c(".senmayo",".SIT")) %>%
  mutate(delta = n.senmayo - n.SIT) %>%
  filter(delta != 0)

# confirm that the total number of cell is the same for both 
full_join(read_tsv("../../out/table/revision/modules_SENESCENCE_shirmer/125_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_senescence_shirmer.tsv") %>%
            filter(signature == "senmayo") %>%
            group_by(sample,cellid,pathology) %>%
            summarise(tot = sum(n)),
          read_tsv("../../out/table/revision/124_enumeration_SIT_WM_CX_harmonySkipIntegAllSoupX_expertAnno_shirmer.tsv") %>%
            group_by(sample,cellid,pathology) %>%
            summarise(tot = sum(n)),by = c("sample","cellid","pathology"),suffix = c(".senmayo",".SIT")) %>%
  mutate(delta = tot.senmayo - tot.SIT) %>%
  filter(delta != 0)


#
df_plot <- full_join(df_senmayo %>%
                       select(sample,pathology,cellid,n,tot,prop),
                     df_SIT %>%
                       select(sample,pathology,cellid,n,tot,prop),by = c("sample","pathology","cellid"),suffix = c(".senmayo",".SIT"))

# it is safe to assume that if one is missing it means it is a 0%
df_plot2 <- df_plot %>%
  mutate(prop.SIT = case_when(is.na(prop.SIT)~0,
                              T~prop.SIT),
         prop.senmayo = case_when(is.na(prop.senmayo)~0,
                                  T~prop.senmayo))


# plot the correlation
# notice that in this case I am discarding the samples for which there is a estimate of 0 senescence
df_plot %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.SIT))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  facet_wrap(~cellid,scales = "free")+theme(strip.background = element_blank())
ggsave(paste0("../../out/image/revision/124_correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_split_shirmer.pdf"),width = 9,height =8)

# notice that in this case I am discarding the samples for which there is a estimate of 0 senescence
df_plot2 %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.SIT))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  facet_wrap(~cellid,scales = "free")+theme(strip.background = element_blank())
ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/124_correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_split_keepZero_shirmer.pdf"),width = 9,height =8)

# plot a globa correlation
df_plot %>%
  filter(!is.na(prop.SIT),
         !is.na(prop.senmayo)) %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.SIT))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/124_correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_global_shirmer.pdf"),width = 5,height =5)

df_plot %>%
  filter(!is.na(prop.SIT),
         !is.na(prop.senmayo)) %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.SIT))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  theme(strip.background = element_blank())+
  scale_x_log10()+
  scale_y_log10()

df_plot %>%
  filter(!is.na(prop.SIT),
         !is.na(prop.senmayo)) %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.SIT))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  theme(strip.background = element_blank())+
  scale_x_sqrt()+
  scale_y_sqrt()

df_plot
cor.test(df_plot$prop.senmayo,df_plot$prop.SIT)

df_plot2 %>%
  filter(!is.na(prop.SIT),
         !is.na(prop.senmayo)) %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.SIT))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/124_correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_global_keepZero_shirmer.pdf"),width = 5,height =5)

df_plot2
cor.test(df_plot2$prop.senmayo,df_plot2$prop.SIT)

df_plot %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.SIT))+
  geom_smooth(method = "lm") +
  geom_point()+theme_bw()+
  facet_wrap(~pathology,scales = "free")+theme(strip.background = element_blank())

# try plot UMAPs for SIT --------------------------------------------------
# pull the annotation derived from the SIT algorithm
df_SIT <- read_tsv("../../out/table/revision/124_meta_SIT_WM_CX_harmonySkipIntegration_AllSoupX_test_shirmer.tsv")

df_SIT %>%
  arrange(SENEQUANTILE) %>%
  mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active"))) %>%
  ggplot(aes(x=curated_1,y=curated_2,col=SENEQUANTILE))+
  geom_point(alpha=0.1,size=0.1)+scale_color_manual(values = c("gray","red"))+
  facet_grid(SENEQUANTILE~pathology)+
  theme_void() +
  guides(color = guide_legend(override.aes = list(size=1,alpha=1)))+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/124_UMAP_SIT_WM_CX_harmonySkipIntegAllSoupX_expertAnno_shirmer.pdf"),width = 7,height = 4)

# balance the number of cells per sample
# what is the minimum nuber of cells per pathology
df_SIT %>% group_by(pathology) %>% summarise(n = n())

set.seed(2144)
df_SIT %>%
  group_by(pathology) %>%
  sample_n(4000) %>%
  arrange(SENEQUANTILE) %>%
  mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active"))) %>%
  ggplot(aes(x=curated_1,y=curated_2,col=SENEQUANTILE))+
  geom_point(alpha=0.1,size=0.1)+scale_color_manual(values = c("gray","red"))+
  facet_grid(SENEQUANTILE~pathology)+
  theme_void() +
  guides(color = guide_legend(override.aes = list(size=1,alpha=1)))+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/124_UMAP_SIT_WM_CX_harmonySkipIntegAllSoupX_expertAnno2_shirmer.pdf"),width = 7,height = 4)
