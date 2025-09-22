# AIM ---------------------------------------------------------------------
# run some preprocessing and QC of the dataset to understand if there are outliers or global trends.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(vsn)
library(hexbin)
library(viridis)
library(pheatmap)
library(PoiClaClu)
library(AnnotationDbi) 
library(AnnotationHub)
library(GGally)
library(RNAseqQC)
library(edgeR)

# read in the data --------------------------------------------------------
ddsHTSeq <- readRDS("../../out/object/201_dds_all.rds")

# remove low expressed genes ----------------------------------------------
# To reduce the noise filter out the lowly expressed genes
# plot the total counts per sample
colSums(counts(ddsHTSeq)) %>%
  data.frame(tot_counts=.) %>%
  rownames_to_column("sample") %>% 
  # add the metadata
  left_join(colData(ddsHTSeq) %>% 
              data.frame(),by = c("sample" = "HUGE.ID")) %>% 
  ggplot(aes(x=ID.subject,y = tot_counts))+geom_col()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = margin(0.5, 0.5, 2, 2, "cm"))
ggsave("../../out/image/202_barplot_tot_count_all.pdf",width = 7,height = 5)

# check the total size of the matrix before filtering fow low expressed genes
nrow(ddsHTSeq)

# remove potential non infirmative genes (lowly expressed genes).
# use the automatic implementation provided by edgeR
ddsHTSeq_filter <- ddsHTSeq[edgeR::filterByExpr(counts(ddsHTSeq), group = colData(ddsHTSeq)$treat),]

# this is the old implementaion 
# ddsHTSeq_filter <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 100, ]

# check the size of the dataset after filtering for lowly expressed genes
nrow(ddsHTSeq_filter)

# save the filtered object
saveRDS(ddsHTSeq_filter,file = "../../out/object/202_dds_all_filter.rds")

# scaling transformation of the data --------------------------------------
# vsd
vds_filter <- vst(ddsHTSeq_filter, blind = T)
head(assay(vds_filter), 3)

# diagnostic plot
meanSdPlot_vsd <- meanSdPlot(assay(vds_filter))
ggsave(plot = meanSdPlot_vsd$gg+theme_bw(),filename = "../../out/image/202_meanSdPlot_vsd_all.pdf",width = 4,height = 4)

# save the object of interest ---------------------------------------------
saveRDS(vds_filter,file = "../../out/object/vds_all_filter.rds")
