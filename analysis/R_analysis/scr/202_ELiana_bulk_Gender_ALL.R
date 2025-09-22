# AIM ---------------------------------------------------------------------
# plot the gender fo the samples

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(UpSetR)
library(gplots)
library(scales)
library(pals)

# plot genes with heatmap -------------------------------------------------
vds_filter <- readRDS(file = "../../out/object/201_dds_all.rds")

gene_id <- rownames(assay(vds_filter)) %in% c("XIST", "USP9Y", "TXLNGY", "EIF1AY", "ZFY", "DDX3Y")
mat <- assay(vds_filter)[gene_id, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat,useNames = TRUE)

# change the rownames to match the symbol 
# rownames(mat2) <- rownames(mat2) %>% 
#   data.frame(ensembl=.) %>% 
#   left_join(DEG_1,by = "ensembl") %>% 
#   pull(symbol) 

# build the annotation object  
meta <- colData(vds_filter) %>% 
  data.frame()

LUT_sample_matrix <- data.frame(sample_matrix = colnames(mat2)) %>%
  left_join(meta,by = c("sample_matrix"="HUGE.ID"))

# sample_ordered <- LUT_sample_matrix$clone

# define a LUT for the clones
# color_id <- alphabet(length(unique(LUT_sample_matrix$clone)))
# check the colors
# show_col(color_id)

# build the named vector
# names(color_id) <- unique(LUT_sample_matrix$clone)

column_ha <- HeatmapAnnotation(gender = LUT_sample_matrix$sex,
                               treat = LUT_sample_matrix$Treatment,
                               group = LUT_sample_matrix$Group,
                               col = list(gender = c("M" = "blue", "F" = "pink"),
                                          group = c("CTRL" = "gray", "MS" = "black"))) 

# change the name of the column in the matrix
colnames(mat2) <- LUT_sample_matrix$ID.subject

# row_ha <- rowAnnotation(class = rep(c("common B_D","offtarget B","offset D"),c(22,3,24)),
#                         col = list(class = c("common B_D" = "violet", "offtarget B" = "yellow","offset D"="brown")))

ht2 <- Heatmap(mat2, 
               name = "exp",
               top_annotation = column_ha, 
               # cluster_rows = F, 
               # col = colorRamp2(c(-2, 0, 1), c("green", "white", "red")),
               # right_annotation = row_ha, 
               # row_split = rep(c(1,2,3,4),c(2,3,4,7)),
               column_title = "gender genes") 

pdf("../../out/image/202_heatmap_GOI_gender_ALL.pdf",width = 8,height = 5) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()