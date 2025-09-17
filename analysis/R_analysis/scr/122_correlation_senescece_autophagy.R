# AIM ---------------------------------------------------------------------
# correlate the autophagy scores with the senescence scores.

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

# load the scores ---------------------------------------------------------
# senescence
folder_sen <-  "../../out/table/revision/modules_SENESCENCE/"
# focus on the autophagy only sigantures
pattern_sig_sen <- paste0(names(readRDS("../../data/signatures/senescence_pathways.rds")),collapse = "|")

file_sen <- dir(folder_sen) %>% 
  str_subset(pattern = "Module") %>%
  str_subset(pattern = pattern_sig_sen)

# they already have the expert annotation suggested by martina
df_modules_sen <- lapply(file_sen, function(x){
  read_tsv(paste0(folder_sen,x))
}) %>% 
  bind_rows() %>% 
  mutate(disease = factor(disease,levels = c("CTRL","MS")))

# autophagy
folder_aut <-  "../../out/table/revision/modules_AUTOPHAGY/"
# focus on the autophagy only sigantures
pattern_sig_aut <- paste0(names(readRDS("../../data/signatures/autophagy_pathways.rds")),collapse = "|")

file_aut <- dir(folder_aut) %>% 
  str_subset(pattern = "Module") %>%
  str_subset(pattern = pattern_sig_aut)

# they already have the expert annotation suggested by martina
df_modules_aut <- lapply(file_aut, function(x){
  read_tsv(paste0(folder_aut,x))
}) %>% 
  bind_rows() %>% 
  mutate(disease = factor(disease,levels = c("CTRL","MS")))


# wrangling ---------------------------------------------------------------
test_sen <- df_modules_sen %>%
  select(barcode,orig.ident,origin,disease,pathology_class,expertAnno.l1,signature,signature_score1) %>%
  filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
  mutate(signature = case_when(signature == "Induces"~"CellAge Induces",
                                signature == "Inhibits"~"CellAge Inhibits",
                                signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
                                signature == "senmayo"~"Senmayo")) %>%
  # average the scores per sample
  group_by(orig.ident,origin,disease,pathology_class,expertAnno.l1,signature) %>%
  summarise(avg_score = mean(signature_score1))
  # pivot_wider(names_from = signature,values_from = avg_score)

# order the table by the test_sen
test_aut <- df_modules_aut %>%
  select(barcode,orig.ident,origin,disease,pathology_class,expertAnno.l1,signature,signature_score1) %>%
  group_by(orig.ident,origin,disease,pathology_class,expertAnno.l1,signature) %>%
  summarise(avg_score = mean(signature_score1))
  # pivot_wider(names_from = signature,values_from = avg_score)

# calculate the correlation matrix
df_full_join <- full_join(test_sen %>%
            pivot_wider(names_from = signature,values_from = avg_score),
          test_aut %>%
            pivot_wider(names_from = signature,values_from = avg_score),
          by = c("orig.ident","origin","disease","pathology_class","expertAnno.l1"),suffix = c(".sen",".aut")) %>%
  ungroup()

# calculate the correlation per score 
# export the matrix for the senescence signatures
mat_cor <- cor(df_full_join %>%
      select(unique(test_sen$signature)),
    df_full_join %>%
      select(unique(test_aut$signature)))

ht <- Heatmap(mat_cor,
        name = "Correlation \nPearson",
        # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
        col = viridis::viridis(option = "turbo",n = 20),
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 8),
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 8),
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE,
        row_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_column_names = T,
        show_row_names = T)

pdf("../../out/image/revision/122_correlation_AutophagySenescence_pearson_matrix.pdf",width = 5,height = 5)
draw(ht,heatmap_legend_side = "left",padding = unit(c(20, 2, 2, 2), "mm"))
dev.off()

# build the dataset for the correlatino plot
df_crossing <- crossing(sig_sen = unique(test_sen$signature),
         sig_aut = unique(test_aut$signature))

# build the scatter plot
df_plot_scatter <- pmap(list(sig_sen = df_crossing$sig_sen,
          sig_aut = df_crossing$sig_aut), function(sig_sen,sig_aut){
            
            df <- df_full_join %>%
              select(orig.ident,origin,disease,pathology_class,expertAnno.l1,sen_score = sig_sen,aut_score = sig_aut) %>%
              mutate(sig_sen_name = sig_sen,
                     sig_aut_name = sig_aut,)
            return(df)
          }) %>%
  bind_rows()

# plot
df_plot_scatter %>%
  ggplot(aes(x=sen_score,y=aut_score))+
  geom_smooth(method = "lm")+
  geom_point(shape=1)+theme_bw()+
  facet_wrap(sig_sen_name~sig_aut_name,scales = "free",ncol=8)+
  theme(strip.background = element_blank())
ggsave("../../out/image/revision/122_correlation_AutophagySenescence_scatter_global.pdf",width = 24,height = 12)

# scatter plot per cell type for senmayo vs regulation of autophagy
df_plot_scatter %>%
  filter(sig_aut_name == "GOBP_REGULATION_OF_AUTOPHAGY") %>%
  filter(sig_sen_name == "Senmayo") %>%
  ggplot(aes(x=sen_score,y=aut_score))+
  geom_smooth(method = "lm")+
  geom_point(shape=1)+theme_bw()+
  facet_wrap(~expertAnno.l1,scales = "free")+
  theme(strip.background = element_blank())
ggsave("../../out/image/revision/122_correlation_AutophagySenescence_scatter_senmayo.pdf",width = 9,height =8)

# save the table
df_plot_scatter %>%
  filter(sig_aut_name == "GOBP_REGULATION_OF_AUTOPHAGY") %>%
  filter(sig_sen_name == "Senmayo") %>%
  write_tsv("../../out/table/revision/122_correlation_AutophagySenescence_scatter_senmayo.tsv")

# use the proportions calculated ------------------------------------------
# join the two datasets pull only the senescence cells from both
df_senmayo <- read_tsv("../../out/table/revision/modules_SENESCENCE/121_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_senescence.tsv") %>%
  filter(signature == "senmayo") %>%
  filter(sen == 1)

read_tsv("../../out/table/revision/modules_SENESCENCE/121_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_senescence.tsv") %>%
  filter(signature == "senmayo") %>%
  filter(sen == 0)

df_autophagy <- read_tsv("../../out/table/revision/modules_AUTOPHAGY/121_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_autophagy.tsv") %>%
  filter(signature == "GOBP_REGULATION_OF_AUTOPHAGY") %>%
  filter(sen == 1)

read_tsv("../../out/table/revision/modules_AUTOPHAGY/121_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_autophagy.tsv") %>%
  filter(signature == "GOBP_REGULATION_OF_AUTOPHAGY") %>%
  filter(sen == 0)

# df_SIT <- read_tsv("../../out/table/revision/124_enumeration_SIT_WM_CX_harmonySkipIntegAllSoupX_expertAnno.tsv") %>%
#   filter(SENEQUANTILE == "YES")

# check the reason for the discrepancies in the two dataset. the call for senescence is not the same in both, therefore it might be absent in one estimate compard to the other

test <- full_join(read_tsv("../../out/table/revision/modules_SENESCENCE/121_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_senescence.tsv") %>%
                    filter(signature == "senmayo") %>%
                    filter(sen == 0) %>%
                    select(signature,orig.ident,expertAnno.l1,pathology_class),
                  read_tsv("../../out/table/revision/modules_AUTOPHAGY/121_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_autophagy.tsv") %>%
                    filter(signature == "GOBP_REGULATION_OF_AUTOPHAGY") %>%
                    filter(sen == 0) %>%
                    select(signature,orig.ident,expertAnno.l1,pathology_class),
                  by = c("orig.ident","expertAnno.l1","pathology_class"),suffix = c(".senmayo",".autophagy"))

test %>%
  filter(is.na(signature.senmayo))

test %>%
  filter(is.na(signature.autophagy))

full_join(df_senmayo %>%
            group_by(orig.ident,origin,pathology_class) %>%
            summarise(n = n()),
          df_autophagy %>%
            group_by(orig.ident,origin,pathology_class) %>%
            summarise(n = n()),by = c("orig.ident","origin","pathology_class"),suffix = c(".senmayo",".autophagy")) %>%
  mutate(delta = n.senmayo - n.autophagy) %>%
  filter(delta != 0)

# confirm that the total number of cell is the same for both 
full_join(read_tsv("../../out/table/revision/modules_SENESCENCE/121_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_senescence.tsv") %>%
            filter(signature == "senmayo") %>%
            group_by(orig.ident,expertAnno.l1,pathology_class) %>%
            summarise(tot = sum(n)),
          read_tsv("../../out/table/revision/modules_AUTOPHAGY/121_enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise_autophagy.tsv") %>%
            filter(signature == "GOBP_REGULATION_OF_AUTOPHAGY") %>%
            group_by(orig.ident,expertAnno.l1,pathology_class) %>%
            summarise(tot = sum(n)),
          by = c("orig.ident","expertAnno.l1","pathology_class"),suffix = c(".senmayo",".autophagy")) %>%
  mutate(delta = tot.senmayo - tot.autophagy) %>%
  filter(delta != 0)

#
df_plot <- full_join(df_senmayo %>%
                       select(orig.ident,origin,pathology_class,expertAnno.l1,n,tot,prop),
                     df_autophagy %>%
                       select(orig.ident,origin,pathology_class,expertAnno.l1,n,tot,prop),by = c("orig.ident","origin","pathology_class","expertAnno.l1"),suffix = c(".senmayo",".autophagy"))

# it is safe to assume that if one is missing it means it is a 0%
df_plot2 <- df_plot %>%
  mutate(prop.autophagy = case_when(is.na(prop.autophagy)~0,
                                    T~prop.autophagy),
         prop.senmayo = case_when(is.na(prop.senmayo)~0,
                                  T~prop.senmayo))


# plot the correlation
# notice that in this case I am discarding the samples for which there is a estimate of 0 senescence
df_plot %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.autophagy))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  facet_wrap(~expertAnno.l1,scales = "free")+theme(strip.background = element_blank())
ggsave(paste0("../../out/image/revision/124_correlation_AutophagySenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_split.pdf"),width = 9,height =8)

# notice that in this case I am discarding the samples for which there is a estimate of 0 senescence
df_plot2 %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.autophagy))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  facet_wrap(~expertAnno.l1,scales = "free")+theme(strip.background = element_blank())
ggsave(paste0("../../out/image/revision/124_correlation_AutophagySenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_split_keepZero.pdf"),width = 9,height =8)

# plot a globa correlation
df_plot %>%
  filter(!is.na(prop.autophagy),
         !is.na(prop.senmayo)) %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.autophagy))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/image/revision/124_correlation_AutophagySenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_global.pdf"),width = 5,height =5)

df_plot %>%
  filter(!is.na(prop.autophagy),
         !is.na(prop.senmayo)) %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.autophagy))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  theme(strip.background = element_blank())+
  scale_x_log10()+
  scale_y_log10()

df_plot %>%
  filter(!is.na(prop.autophagy),
         !is.na(prop.senmayo)) %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.autophagy))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  theme(strip.background = element_blank())+
  scale_x_sqrt()+
  scale_y_sqrt()

df_plot
cor.test(df_plot$prop.senmayo,df_plot$prop.autophagy)

df_plot2 %>%
  filter(!is.na(prop.autophagy),
         !is.na(prop.senmayo)) %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.autophagy))+
  geom_smooth(method = "lm") +
  geom_point(alpha=0.5)+theme_bw()+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/image/revision/124_correlation_AutophagySenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_global_keepZero.pdf"),width = 5,height =5)

df_plot2
cor.test(df_plot2$prop.senmayo,df_plot2$prop.autophagy)

df_plot %>%
  ggplot(aes(x=prop.senmayo,
             y=prop.autophagy))+
  geom_smooth(method = "lm") +
  geom_point()+theme_bw()+
  facet_wrap(~pathology_class,scales = "free")+theme(strip.background = element_blank())

# # correlation matrix ------------------------------------------------------
# test_correlation <- cor(exp_avg_wide %>%
#                           column_to_rownames("group_id"))
# # plot as heatmap
# pdf("../../out/plot/correlation_stemness_pearson_matrix.pdf",width = 5,height = 4)
# Heatmap(test_correlation,
#         name = "Correlation \nPearson",
#         col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
#         row_names_side = "right",
#         row_names_gp = gpar(fontsize = 8),
#         column_names_side = "bottom",
#         column_names_gp = gpar(fontsize = 8),
#         row_dend_reorder = FALSE,
#         column_dend_reorder = FALSE,
#         row_title_gp = gpar(fontsize = 10, fontface = "bold"),
#         column_title_gp = gpar(fontsize = 10, fontface = "bold"),
#         show_column_names = F,
#         show_row_names = T)
# dev.off()