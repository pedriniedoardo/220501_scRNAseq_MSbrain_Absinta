# AIM ---------------------------------------------------------------------
# compare the subclustering of AST with the original object

# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
library(harmony)
library(ggrepel)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
# read in the cells from run 02
# remove the cells from the TBHP samples
sobj_tot <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
DimPlot(sobj_tot,raster = T,group.by = "expertAnno.l1")

# read in the subcluster
sobj_AST <- readRDS("../../out/object/131_AST_subcluster_HarmonySample.rds")
DimPlot(sobj_AST,group.by = "RNA_snn_res.0.1",raster = T)

# wrangling ---------------------------------------------------------------
# add the rename cluster Id from the res 0.1 to the original object
meta_test_subcluster <- sobj_AST@meta.data %>%
  rownames_to_column("barcode") %>%
  dplyr::select(clu_AST_sub = RNA_snn_res.0.1,barcode)

# add the annotation from the cubsluster to the original object
sobj_tot@meta.data$test_cluster <- sobj_tot@meta.data %>%
  rownames_to_column("barcode") %>%
  left_join(meta_test_subcluster,by = c("barcode")) %>%
  mutate(test_cluster = case_when(is.na(clu_AST_sub) ~ expertAnno.l1,
                                  # is.na(clu_MG_sub) & expertAnno.l1 != "MG"~expertAnno.l1,
                                  # is.na(clu_MG_sub) & expertAnno.l1 == "MG"~"MG_cleaned",
                                  T~clu_AST_sub)) %>%
  pull(test_cluster)

# make the plotting
shortlist_features_list_long <- list(
  IMMUNE = c("LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC","AIF1","HLA-DRA","TYROBP"),
  B_CELLS = c("IGHG1", "CD38"),
  T_CELLS =  c("SKAP1", "CD8A", "CD2"),
  OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
  NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
  PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"),
  SCHWANN = c("PMP22","MPZ","PRX"),
  EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
  STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
# source image
test_long_run03 <- DotPlot(sobj_tot,
                           features = shortlist_features_list_long,
                           dot.scale = 8,
                           cluster.idents = T,
                           group.by = "expertAnno.l1") +
  RotatedAxis() +
  labs(title = "full WM+CX expertAnno.l1")+
  theme(strip.text = element_text(angle = 90))

# generate the ggplot object
# force the order suggested by matrina
df_test_run03 <- lapply(shortlist_features_list_long,function(x){
  test_long_run03$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

# plot mapping radius
test_long_run03_alt <- df_test_run03 %>%
  # force the order
  mutate(id = factor(id,levels = c("IMM","LYM","INH NEU","EXC NEU","AST","OLIGO","VAS","OPC","EPENDYMA"))) %>% 
  mutate(cell_type = factor(cell_type,levels = c("IMMUNE","B_CELLS","T_CELLS","NEURONS","ASTRO","OLIGOLINEAGE","ENDO","PERICYTE","SCHWANN","EPENDYMA","STROMAL"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90),strip.text = element_text(hjust = 0,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")

#
test_long_subAST <- DotPlot(sobj_AST,
                            features = shortlist_features_list_long,
                            dot.scale = 8,
                            cluster.idents = T,
                            group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "AST subset res 0.1")+
  theme(strip.text = element_text(angle = 90))

# generate the ggplot object
# force the order suggested by matrina
df_test_subAST <- lapply(shortlist_features_list_long,function(x){
  test_long_subAST$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

# plot mapping radius
test_long_subMG_alt <- df_test_subAST %>%
  # force the order
  # mutate(id = factor(id,levels = c("IMM","LYM","INH NEU","EXC NEU","AST","OLIGO","OPC","AST","EPENDYMA"))) %>% 
  mutate(cell_type = factor(cell_type,levels = c("IMMUNE","B_CELLS","T_CELLS","NEURONS","ASTRO","OLIGOLINEAGE","ENDO","PERICYTE","SCHWANN","EPENDYMA","STROMAL"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90),strip.text = element_text(hjust = 0,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")

#
test_long_run032 <- DotPlot(sobj_tot,
                            features = shortlist_features_list_long,
                            dot.scale = 8,
                            cluster.idents = T,
                            group.by = "test_cluster") +
  RotatedAxis() +
  labs(title = "WM+CX with AST subcluster")+
  theme(strip.text = element_text(angle = 90))

# generate the ggplot object
# force the order suggested by matrina
df_test_run032 <- lapply(shortlist_features_list_long,function(x){
  test_long_run032$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

# plot mapping radius
test_long_run032_alt <- df_test_run032 %>%
  # force the order
  mutate(id = factor(id,levels = c(0:4,"IMM","LYM","INH NEU","EXC NEU","VAS","OLIGO","OPC","EPENDYMA"))) %>%
  mutate(cell_type = factor(cell_type,levels = c("IMMUNE","B_CELLS","T_CELLS","NEURONS","ASTRO","OLIGOLINEAGE","ENDO","PERICYTE","SCHWANN","EPENDYMA","STROMAL"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90),strip.text = element_text(hjust = 0,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")

ggsave(plot = (test_long_run03_alt / test_long_subMG_alt / test_long_run032_alt),"../../out/image/131_DotplotLong_AST_subcluster_test.pdf",width = 20,height = 18)
