# AIM ---------------------------------------------------------------------
# run some preprocessing and QC of the dataset to understand if there are outliers or global trends.
# update data analysis
# this version of the analysis uses the trt vs naive data

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
library(pals)
library(scales)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
ddsHTSeq <- readRDS("../../out/object/201_dds_trt_update.rds")

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
ggsave("../../out/image/202_barplot_trt_count_trt_update.pdf",width = 7,height = 5)

# check the total size of the matrix before filtering fow low expressed genes
nrow(ddsHTSeq)

# remove potential non infirmative genes (lowly expressed genes).
# use the automatic implementation provided by edgeR
ddsHTSeq_filter <- ddsHTSeq[edgeR::filterByExpr(counts(ddsHTSeq), group = colData(ddsHTSeq)$Group),]

# this is the old implementaion 
# ddsHTSeq_filter <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 100, ]

# check the size of the dataset after filtering for lowly expressed genes
nrow(ddsHTSeq_filter)

# save the filtered object
saveRDS(ddsHTSeq_filter,file = "../../out/object/202_dds_trt_filter_update.rds")

# scaling transformation of the data --------------------------------------
# vsd
# save also the value for the unfiltered values
vds <- vst(ddsHTSeq, blind = T)
vds_filter <- vst(ddsHTSeq_filter, blind = T)

head(assay(vds), 3)
head(assay(vds_filter), 3)

# diagnostic plot
meanSdPlot_vsd <- meanSdPlot(assay(vds_filter))
ggsave(plot = meanSdPlot_vsd$gg+theme_bw(),filename = "../../out/image/202_meanSdPlot_vsd_trt_update.pdf",width = 4,height = 4)

# save the object of interest ---------------------------------------------
saveRDS(vds,file = "../../out/object/vds_trt_unfilter_update.rds")
saveRDS(vds_filter,file = "../../out/object/vds_trt_filter_update.rds")


# plot cluster alternative ------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/image/202_heatmap_trt_clustering.pdf",width = 7,height = 5)
set.seed(2)
plot_sample_clustering(vds_filter,
                       anno_vars = c("sex","Treatment","Group"),
                       distance = "euclidean")
dev.off()

# PCA plot ----------------------------------------------------------------
plot_vsd <- plotPCA(vds_filter,
                    intgroup = c("sex","Treatment","Group")) +
  theme_bw()

plot_vsd$data %>%
  ggplot(aes(x=PC1,y=PC2,label = name)) +
  # ggplot(aes(x=PC1,y=PC2,col=BMP_treat,shape=condition)) +
  geom_point(aes(col=Treatment,shape=sex),size =3) +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave("../../out/image/202_PCA_vsd_trt_colorGroup.pdf",width = 10,height = 8)

# pull more PC
rv <- rowVars(assay(vds_filter),useNames = TRUE)
# select the ntop genes by variance
select_var <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
test <- prcomp(t(assay(vds_filter)[select_var,]))$x %>% 
  data.frame() %>% 
  rownames_to_column("sample")

# plot more PC by condition
left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample")) %>%
  # ggpairs(columns = 5:14,ggplot2::aes(colour=condition),upper = "blank")+
  ggpairs(columns = 6:11,ggplot2::aes(colour=Treatment),upper = "blank") +
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/202_PCA_panel_vsd_trt_treat.pdf",width = 10,height = 10)

# explore pc score by metadata fo the samples
test2 <- left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample"))

test_df1 <- test2 |> 
  dplyr::select(sex,Treatment,name) |> 
  # dplyr::rename(approx.time = Approx.Time.between.collection.and.processing) |> 
  pivot_longer(names_to = "var_1",values_to = "value_1",-c(name))

test_df2 <- test2 |> 
  dplyr::select(name,PC1:PC6) |>
  pivot_longer(names_to = "var_2",values_to = "value_2",-c(name))

# martina asked to label the sample RNA_MRA_04
left_join(test_df1,test_df2,by=c("name")) |> 
  mutate(comparison = paste0(var_1,"_vs_",var_2)) |> 
  mutate(test = case_when(name %in% c("RNA_MRA_04")~"RNA_MRA_04",
                          T ~ "other")) |>
  ggplot(aes(x=value_1,y=value_2)) +
  facet_wrap(~comparison,scales = "free",ncol=6) +
  # facet_grid(var_1~var_2,scales = "free_x",drop = T) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2),aes(col=test))+
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/202_panel_trt_metadata_PC.pdf",width = 15,height = 6)

# plot gender -------------------------------------------------------------
GOI <- c(
  "XIST",   # Expressed in XX individuals to inactivate one X chromosome
  # "WNT4",   # Actively promotes ovary development
  # "RSPO1",  # Works with WNT4 to promote ovary development
  # "FOXL2",   # Maintains ovary function and suppresses the male pathway
  # "SRY",      # Master switch for male sex determination
  "DDX3Y",    # Sperm production
  "USP9Y",    # Sperm production
  "EIF1AY"   # Sperm production
  # "NLGN4Y",   # Involved in brain synapse formation
  # "AMELY",     # Involved in tooth enamel formation (used in forensics)
  # "PCDH11Y"  # Involved in brain cell recognition
  
)
# gene_id <- rownames(assay(vds)) %in% df_GOI$GENEID
# mat <- assay(vds)[gene_id, ]
mat <- assay(vds)[GOI, ]
mat2 <- (mat - rowMeans(mat))/rowSds(mat,useNames = TRUE)

# add the gene symbol to the matrix
# rownames(mat2) <- df_GOI$SYMBOL

# build the annotation object  
meta <- colData(vds) %>% 
  data.frame()

LUT_sample_matrix <- data.frame(sample_matrix = colnames(mat2)) %>%
  left_join(meta,by = c("sample_matrix"="HUGE.ID"))

# sample_ordered <- LUT_sample_matrix$clone

# define a LUT for the clones
# color_id <- alphabet(length(unique(LUT_sample_matrix$Group)))
# check the colors
# show_col(color_id)

# build the named vector
# names(color_id) <- unique(LUT_sample_matrix$clone)

column_ha <- HeatmapAnnotation(# clone = LUT_sample_matrix$clone,
  gender = LUT_sample_matrix$sex,
  treat = LUT_sample_matrix$Treatment,
  group = LUT_sample_matrix$Group,
  col = list(# clone = color_id,
    gender = c("M" = "blue", "F" = "pink"))) 

# change the name of the column in the matrix
# colnames(mat2) <- LUT_sample_matrix$test

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

pdf("../../out/image/202_heatmap_GOI_trt_gender.pdf",width = 8,height = 5) 
draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
dev.off()
