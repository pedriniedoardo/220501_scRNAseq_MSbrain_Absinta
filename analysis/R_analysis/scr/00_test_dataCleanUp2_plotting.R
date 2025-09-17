# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggraph)
library(igraph)
library(data.tree)
library(HGNChelper)

# read in the data --------------------------------------------------------
# result of the cleanup suggested by martina
scobj <- readRDS("../../out/object/00_test_dataCleanUP_test.rds")
DimPlot(scobj,raster = T,reduction = "curated",label = T,split.by = "origin_alt")
DimPlot(scobj,raster = T,reduction = "curated",label = T)

# save the panel martina uses for the annotation
shortlist_features_list_long <- list(
  IMMUNE = c("IGHG1","CD38","CD8A","CD2","SKAP1","LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
  NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
  PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long2 <- DotPlot(scobj, features = shortlist_features_list_long, dot.scale = 8,cluster.idents = T)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# theme(axis.text.x = element_text(angle = 45,hjust = 1))
test_long2
ggsave("../../out/image/ManualClean/Dotplot_test_dataCleanUP_test.pdf",width = 23,height = 4)

# -------------------------------------------------------------------------
# load the original dataset and compare the cluster identity
# result of the cleanup suggested by martina
scobj2 <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegration_AllSoupX_test.rds")
DimPlot(scobj2,raster = T,label = T)

# are the clusters identities the same?
test <- left_join(x = 
            scobj@meta.data %>%
            rownames_to_column("barcode") %>%
            dplyr::select(barcode,seurat_clusters),
          y = scobj2@meta.data %>%
            rownames_to_column("barcode") %>%
            dplyr::select(barcode,seurat_clusters),suffix = c(".ref",".martina"),by="barcode")
  
table(test$seurat_clusters.ref,test$seurat_clusters.martina)
