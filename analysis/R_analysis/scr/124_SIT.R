##########################################################
############## SENESCENCE INDEX TOOL (SIT) ###############
##########################################################

# libraries ---------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(homologene)
library(RColorBrewer)
library(scales)
library(msigdbr)
library(clusterProfiler)
library(tidyverse)
library(GSVA)
library(limma)
library(cowplot)

# library(org.Mm.eg.db)

# read in the data --------------------------------------------------------
Seurat.object <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
DimPlot(Seurat.object,raster = T)

# run the sample scoring approach -----------------------------------------

########## SENESCENCE SIGNATURE ######
# reference markers from the authors
senescence_marker_mouse <- c("Cdkn2a","Cdkn1a",  "Serpine1" ,"Cdkn1b", "Cdkn2d", "Cdkn2b")

# convert the human gene neames
df2 <- homologene(senescence_marker_mouse, inTax = 10090, outTax = 9606) %>% 
  setNames(c("mouse_gene","human_gene","mouse_ID","human_ID"))
df2

senescence_marker <- df2$human_gene

# VlnPlot(Seurat.object, senescence_marker,raster = T)

# add the score fo the subset of markers
Seurat.object <- AddModuleScore(Seurat.object, features = list(senescence_marker),name = "Senescence_signature")
head(Seurat.object@meta.data)

# library(scales)
FeaturePlot(Seurat.object, features = "Senescence_signature1",raster = T,order = T)+
  scale_color_gradientn(colours =  rev(brewer.pal(n = 11, name = "RdBu")), 
                        values = rescale(c(min(Seurat.object$Senescence_signature1),0,max(Seurat.object$Senescence_signature1))),guide = "colorbar", limits = c(min(Seurat.object$Senescence_signature1), max(Seurat.object$Senescence_signature1)))

########## CELL CYCLE SIGNATURE ######
# get the specific pathways
# get a shortlist of the terms
m_df <- msigdbr(species = "Homo sapiens",category = "C2")
head(m_df)

# filter only the terms fo interest and pull the genes
GSred <- m_df %>% 
  filter(gs_name %in% c("KEGG_CELL_CYCLE","REACTOME_CELL_CYCLE" ,"WP_CELL_CYCLE")) %>% 
  split(f = .$gs_name) %>% 
  map(function(x){
    x %>% 
      pull(gene_symbol)
  })

str(GSred)

# run the enrichment analysis per cell
# library(GSVA)

# -------------------------------------------------------------------------
# this is too mem intensive
# exprdf <- as.matrix(Seurat.object@assays$RNA@data)
# gsvascore1 <- GSVA::gsva(expr = exprdf, GSred[1],  mx.diff=T, method="ssgsea", parallel.sz=24)

# test use less cores to reduce RAM memory
gsvascore <- GSVA::gsva(expr = Seurat.object@assays$RNA@data, GSred,  mx.diff=T, method="ssgsea", parallel.sz=4)
# saveRDS(gsvascore,"../../out/object/revision/124_test_gsvascore.rds")
# free the memory
gc()
gsvascore <- normalizeQuantiles(gsvascore)
# save the object
saveRDS(gsvascore,"../../out/object/revision/124_gsvascore.rds")
# gsvascore <- readRDS("../../out/object/gsvascore.rds")

str(gsvascore)

colsnew <- c((brewer.pal(n = 9, name = "PiYG")))
colsnew[5] <- "#f7f4b0"

Seurat.object$KEGG_CELL_CYCLE <- gsvascore["KEGG_CELL_CYCLE", colnames(Seurat.object)]
Seurat.object$Cell_cycle_arrest_signatureK <- (-(Seurat.object$KEGG_CELL_CYCLE))
FeaturePlot(Seurat.object, features = "Cell_cycle_arrest_signatureK",raster = T)+
  scale_color_gradientn(colours =  colsnew,
                        guide = "colorbar",
                        limits = c(min(Seurat.object$Cell_cycle_arrest_signatureK), max(Seurat.object$Cell_cycle_arrest_signatureK)))

# Reactome signature add the scores to the object
Seurat.object$REACTOME_CELL_CYCLE <- gsvascore["REACTOME_CELL_CYCLE", colnames(Seurat.object)]
Seurat.object$Cell_cycle_arrest_signatureR <- -(Seurat.object$REACTOME_CELL_CYCLE)
FeaturePlot(Seurat.object, features = "Cell_cycle_arrest_signatureR")+
  scale_color_gradientn(colours =  colsnew,
                        guide = "colorbar",
                        limits = c(min(Seurat.object$Cell_cycle_arrest_signatureR), max(Seurat.object$Cell_cycle_arrest_signatureR)))

# WP signature add the scores to the object
Seurat.object$WP_CELL_CYCLE <- gsvascore["WP_CELL_CYCLE", colnames(Seurat.object)]
Seurat.object$Cell_cycle_arrest_signatureWP <- -(Seurat.object$WP_CELL_CYCLE)
FeaturePlot(Seurat.object, features = "Cell_cycle_arrest_signatureWP")+
  scale_color_gradientn(colours =  colsnew,
                        guide = "colorbar",
                        limits = c(min(Seurat.object$Cell_cycle_arrest_signatureWP), max(Seurat.object$Cell_cycle_arrest_signatureWP)))

########  SENESCENCE SCORE ########
senenorm <- scales::rescale(x=Seurat.object$Senescence_signature1, c(0,1))
keggnorm <- scales::rescale(x=Seurat.object$Cell_cycle_arrest_signatureK, c(0,1))
reactomenorm <- scales::rescale(x=Seurat.object$Cell_cycle_arrest_signatureR, c(0,1))
wpnorm <- scales::rescale(x=Seurat.object$Cell_cycle_arrest_signatureWP, c(0,1))

colors <- rainbow(16, start=0.1, end=0.9)

Seurat.object$Senescence_scoreK <- senenorm+keggnorm
FeaturePlot(Seurat.object, features =  "Senescence_scoreK",raster = T,order=T)+ 
  scale_color_gradientn(colours =  colors)

Seurat.object$Senescence_scoreWP <- senenorm+wpnorm
FeaturePlot(Seurat.object, features =  "Senescence_scoreWP",raster = T,order=T)+ 
  scale_color_gradientn(colours =  colors)

Seurat.object$Senescence_scoreR <- senenorm+reactomenorm
FeaturePlot(Seurat.object, features =  "Senescence_scoreR",raster = T,order=T)+ 
  scale_color_gradientn(colours =  colors)

##########  Identification of senescent cells
quantile <- cbind(as.character(cut(as.numeric(as.character(Seurat.object$Senescence_scoreK)),as.numeric(quantile(Seurat.object$Senescence_scoreK)))),
                  as.character(cut(as.numeric(as.character(Seurat.object$Senescence_scoreR)),as.numeric(quantile(Seurat.object$Senescence_scoreR)))),
                  as.character(cut(as.numeric(as.character(Seurat.object$Senescence_scoreWP)),as.numeric(quantile(Seurat.object$Senescence_scoreWP)))))

# rename the levels of the quantiles using the level nomenclature
df_quantiles <- lapply(1:ncol(quantile),function(x){
  data.frame(x = factor(quantile[,x],labels = c("quantile1","quantile2","quantile3","quantile4")))
}) %>% 
  bind_cols() %>% 
  # rename the columns
  setNames(c("Senescence_scoreKQUANTILE", "Senescence_scoreRQUANTILE", "Senescence_scoreWPQUANTILE")) %>% 
  # add the quantiles
  mutate(barcode = colnames(Seurat.object))

head(df_quantiles)

# ref metadata
df_meta_ref <- Seurat.object@meta.data

# merge the metadata
df_meta_full <- df_meta_ref %>% 
  rownames_to_column("barcode") %>% 
  # add the quantile information
  left_join(df_quantiles,by = "barcode") %>% 
  # add the criterion of selection for the senescent cells. all the cells should be in the top quantiles for all the signatures
  mutate(SENEQUANTILE = case_when(Senescence_scoreKQUANTILE == "quantile4" &
                                    Senescence_scoreRQUANTILE == "quantile4" &
                                    Senescence_scoreWPQUANTILE == "quantile4" ~ "YES",
                                  T ~ "NO"))

# check
df_meta_full %>%
  dplyr::select(barcode,Senescence_scoreK,Senescence_scoreWP,Senescence_scoreR,Senescence_scoreKQUANTILE,Senescence_scoreRQUANTILE,Senescence_scoreWPQUANTILE) %>% 
  pivot_longer(names_to = "signature",values_to = "score",c(Senescence_scoreK,Senescence_scoreWP,Senescence_scoreR)) %>% 
  pivot_longer(names_to = "quantile_signature",values_to = "quantile",c(Senescence_scoreKQUANTILE,Senescence_scoreRQUANTILE,Senescence_scoreWPQUANTILE)) %>% 
  # match the quantile_signature and the signature
  mutate(quantile_signature = str_sub(quantile_signature,start = 1,end = -9)) %>% 
  filter(signature == quantile_signature) %>% 
  ggplot(aes(x=factor(1),y=score))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(col=quantile),position=position_jitter(width = 0.2),alpha=0.5)+facet_wrap(~signature)+theme_bw()+theme(strip.background = element_blank())

# update the metadata in the seurat object
Seurat.object@meta.data <- df_meta_full %>% 
  column_to_rownames("barcode")
# plotting
DimPlot(Seurat.object,raster = T,order = T, group.by = "SENEQUANTILE",cols= c("NO"="#F15A2B" , "YES"="#1D75BC"))+ ggtitle("SENESCENCE_QUANTIFICATION: quantile")

DimPlot(Seurat.object,raster = T, group.by = "SENEQUANTILE",cols= c("NO"="#F15A2B" , "YES"="#1D75BC"),split.by = "pathology_class")+ ggtitle("SENESCENCE_QUANTIFICATION: quantile")

test_NO <- subset(Seurat.object,subset = SENEQUANTILE == "NO")
test_YES <- subset(Seurat.object,subset = SENEQUANTILE == "YES")

DimPlot(test_YES,raster = T, group.by = "SENEQUANTILE",cols= c("NO"="#F15A2B" , "YES"="#1D75BC"),split.by = "pathology_class")+ ggtitle("SENESCENCE_QUANTIFICATION: quantile")

# extract the metadata to be more flexible with the plotting
# save the current meta add also the coordinates of the UMAP
df_umap <- Seurat.object@reductions$umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column()

df_meta <- Seurat.object@meta.data %>% 
  data.frame() %>% 
  rownames_to_column()

df_meta_full <- left_join(df_umap,df_meta,"rowname")

# full map
df_meta_full %>%
  arrange(SENEQUANTILE) %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=SENEQUANTILE))+geom_point(size = 0.1,alpha=0.5)+theme_cowplot()+
  scale_color_manual(values = c("NO"="#F15A2B" , "YES"="#1D75BC"))+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))
ggsave("../../out/image/revision/124_data.combined_WM_CX_harmonySkipIntegration_AllSoupX_SIT.pdf",height = 10,width = 11)
ggsave("../../out/image/revision/124_data.combined_WM_CX_harmonySkipIntegration_AllSoupX_SIT.png",height = 10,width = 11,bg="white")

# split the map by both pathology condition and being senescence or not
df_meta_full %>%
  arrange(SENEQUANTILE) %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=SENEQUANTILE))+geom_point(size = 0.1,alpha=0.5)+theme_cowplot()+
  scale_color_manual(values = c("NO"="#F15A2B" , "YES"="#1D75BC"))+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1))) +
  facet_grid(SENEQUANTILE~pathology_class) + theme(strip.background = element_blank())
ggsave("../../out/image/revision/124_data.combined_WM_CX_harmonySkipIntegration_AllSoupX_SIT_split.pdf",height = 10,width = 41)
ggsave("../../out/image/revision/124_data.combined_WM_CX_harmonySkipIntegration_AllSoupX_SIT_split.png",height = 10,width = 41,bg="white")

# count the senescent cells per sample
df_summary_sample <- df_meta_full %>% 
  group_by(orig.ident,origin,pathology_class,SENEQUANTILE) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(orig.ident,origin,pathology_class) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop_senescence = n/tot)
write_tsv(df_summary_sample,"../../out/table/revision/124_summary_SIT_sample.tsv")

df_summary_sample %>% 
  ggplot(aes(x=pathology_class,y=prop_senescence))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  facet_wrap(~SENEQUANTILE,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/revision/124_data.combined_WM_CX_harmonySkipIntegration_AllSoupX_SIT_propSample.pdf",height = 4,width = 8)
ggsave("../../out/image/revision/124_data.combined_WM_CX_harmonySkipIntegration_AllSoupX_SIT_propSample.png",height = 4,width = 8,bg="white")

# do the same per cluster
df_summary_cluster <- df_meta_full %>% 
  group_by(orig.ident,RNA_snn_res.0.1,origin,pathology_class,SENEQUANTILE) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(orig.ident,RNA_snn_res.0.1,origin,pathology_class) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop_senescence = n/tot)
write_tsv(df_summary_cluster,"../../out/table/revision/124_summary_SIT_cluster.tsv")

df_summary_cluster %>%
  filter(SENEQUANTILE=="YES") %>% 
  ggplot(aes(x=pathology_class,y=prop_senescence))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  facet_wrap(~RNA_snn_res.0.1,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/revision/124_data.combined_WM_CX_harmonySkipIntegration_AllSoupX_SIT_propCluster.pdf",height = 12,width = 12)
ggsave("../../out/image/revision/124_data.combined_WM_CX_harmonySkipIntegration_AllSoupX_SIT_propCluster.png",height = 12,width = 12,bg="white")

# use the cell type imputation from SCType
# df_meta_full <- df_meta_full %>%
#   mutate(cell_id = case_when(seurat_clusters %in% c(0,4,14)~"OLIGO",
#                              seurat_clusters %in% c(9)~"OPC",
#                              seurat_clusters %in% c(3,12)~"ASTRO",
#                              seurat_clusters %in% c(5)~"IMMUNE",
#                              seurat_clusters %in% c(13)~"LYM",
#                              seurat_clusters %in% c(11)~"VAS",
#                              seurat_clusters %in% c(1, 2, 10,6)~"EXC NEU",
#                              seurat_clusters %in% c(7,8)~"INH NEU"))

df_summary_cellID <- df_meta_full %>% 
  group_by(orig.ident,expertAnno.l1,origin,pathology_class,SENEQUANTILE) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(orig.ident,expertAnno.l1,origin,pathology_class) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop_senescence = n/tot)
write_tsv(df_summary_cellID,"../../out/table/revision/124_summary_SIT_cellID.tsv")

df_summary_cellID %>%
  filter(SENEQUANTILE=="YES") %>% 
  ggplot(aes(x=pathology_class,y=prop_senescence))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  facet_wrap(~expertAnno.l1,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/revision/124_data.combined_WM_CX_harmonySkipIntegration_AllSoupX_SIT_propCellID.pdf",height = 9,width = 9)
ggsave("../../out/image/revision/124_data.combined_WM_CX_harmonySkipIntegration_AllSoupX_SIT_propCellID.png",height = 9,width = 9,bg="white")

# VlnPlot(Seurat.object, features = c("Senescence_signature", "Senescence_scoreR", "Senescence_scoreK", "Senescence_scoreWP"), 
#         ncol = 2, group.by = "SENEQUANTILE",cols= c("NO"="#F15A2B" , "YES"="#1D75BC"))

# save the table of the meta after SIT scoring
df_meta_full %>%
  write_tsv("../../out/table/revision/124_meta_SIT_WM_CX_harmonySkipIntegration_AllSoupX_test.tsv")
