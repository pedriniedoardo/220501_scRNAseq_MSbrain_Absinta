# AIM ---------------------------------------------------------------------
# plot some genes of interst
# in particular this is to confirm the trend of the common genes

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(UpSetR)
library(gplots)

# plot genes dotplot ------------------------------------------------------
# save the filtered object
# read_tsv("../../out/table/DE_treatvsCSFctrl_pseudobulk_MG_shr.tsv") %>%
#   filter(conditionVsCSFctrl == "CSF.MS.24h_shr") %>%
#   filter(symbol %in% c("CSF1R","CTSA","CTSB","TREM2"))
data <- readRDS("../../out/object/201_dds_all_update.rds") %>%
  DESeq()
# data <- readRDS("../../out/object/202_dds_all_filter_update.rds") %>%
#   DESeq()

# pull the list of degs reported by eliana
test_DEGS <- readxl::read_xlsx("../../data/signatures/eliana_endo_ms_update/DEGs_MS_ctrl.xlsx")

# rank them by top absolute logFC and pick the topmost disregulated
test_DEGS %>%
  arrange(desc(logCPM)) %>% print(n=40)
  arrange(desc(abs(logFC)))

test_DEGS %>%
  dplyr::filter(gene %in% c("ABCA4"))

# correct the HUGE-ID nameing
lut <- colData(data) %>%
  as.data.frame()
  # mutate(sample = str_replace_all(HUGE.ID,pattern = "-","."))

# GOI <- c(sig_B$symbol,sig_D$symbol)
GOI <- c("SERPINA1","PROX1","PODXL","IGFBP4")

# GOI <- subset_genes
# are all the GOI in the table
# sum(!GOI %in% rownames(test))
sum(!GOI %in% rownames(data))

MR <- counts(data,normalized=F)%>%
  as.data.frame() %>%
  rownames_to_column("symbol") %>%
  pivot_longer(names_to = "sample",values_to = "exp",-symbol)%>%
  group_by(sample)%>%
  summarise(MR = sum(exp)/10^6)

# plot the data following the methods implemented in the plotCounts funciton from DESeq2
# Normalized counts plus a pseudocount of 0.5 are shown by default.
counts(data,normalized=T) %>%
  as.data.frame() %>%
  rownames_to_column("symbol") %>%
  # filter(symbol %in% GOI) %>%
  filter(symbol %in% "ABCA4") %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = c("sample" = "HUGE.ID")) %>%
  mutate(count_norm_adj = count + 0.5)%>%
  ggplot(aes(x=Group,y = count_norm_adj))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6)+facet_wrap(~symbol,scales = "free")+scale_y_log10()+ theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("../../out/image/boxplot_GOI.pdf",width = 6,height = 6)

counts(data,normalized=T)%>%
  as.data.frame()%>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = c("sample" = "HUGE.ID")) %>%
  mutate(count_norm_adj = count + 0.5) %>%
  # ggplot(aes(x=Group,y = count_norm_adj,label=clone))+
  ggplot(aes(x=Group,y = count_norm_adj))+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6)+
  facet_wrap(~symbol,scales = "free") + scale_y_log10() + theme_bw() +
  # geom_text_repel()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("../../out/image/boxplot_GOI_label.pdf",width = 6,height = 6)

# plot full panel of marker genes -----------------------------------------
# plot the full panel of marker genes.
# load the full vst transformed values of expression. use the indiltered values
vst <- readRDS(file = "../../out/object/vds_all_unfilter_update.rds")

# pull the matrix of expresison
mat_filter <- assay(vst) %>%
  data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

setdiff(test_DEGS$gene,rownames(vst))

GOI <- test_DEGS$gene[!test_DEGS$gene %in% c("UTS2","KRTAP2-3","IL2RA")]

# subset the DEGs only
mat <- mat_filter[GOI, ]

# z-scale the matrix
mat2 <- (mat - rowMeans(mat))/rowSds(mat)

# change the rownames with gene symbols
# rownames(mat2) <- df_DEG_1$SYMBOL
#

# confirm the scalinng per gene
# mat2 %>%
#   as.data.frame() %>%
#   rownames_to_column("gene") %>%
#   pivot_longer(names_to = "sample",values_to = "zscore",-gene) %>%
#   group_by(gene) %>%
#   summarise(avg = mean(zscore),
#             sd = sd(zscore)) %>%
#   arrange(desc(avg))

sample_ordered <- data.frame(HUGE.ID = colnames(mat2)) %>%
  left_join(vst@colData %>%
              data.frame() %>%
              mutate(HUGE.ID = str_replace_all(HUGE.ID,pattern = "-",".")),by="HUGE.ID")


# update the column name of the matrix
colnames(mat2) <- sample_ordered$ID.subject

column_ha <- HeatmapAnnotation(disease = sample_ordered$Group,
                               treat = sample_ordered$Treatment,
                               gender = sample_ordered$sex,
                               col = list(treat = c("ctrl" = "red", "Naive"="orange","TRT" = "cyan"),
                                          disease = c("CTRL" = "gray", "MS"="black"),
                                          gender = c("M" = "blue", "F" = "pink"))) 

ht2 <- Heatmap(mat2, 
               name = "exp", 
               column_title = "treat",
               row_names_gp = gpar(fontsize = 5),
               top_annotation = column_ha, show_row_names = F
               # cluster_rows = F, 
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7))
) 
# pdf("../../out/plot/heatmap_DEG_cytokine.pdf",width = 5,height = 6) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
# dev.off()