# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)

# read the data -----------------------------------------------------------
# read in the dataset
data.combined <- readRDS("../../data/schirmer.rds")
# make sure the RNA assay is loaded
DefaultAssay(data.combined) <- "RNA"

# add more annotations
# add the macro annotation from the cellid, based on martina's paper
data.combined$cellid <- data.combined@meta.data %>%
  mutate(cellid = case_when(seurat_clusters %in% c(0,4,10,18)~"OLIGO",
                            seurat_clusters %in% c(14)~"VAS",
                            seurat_clusters %in% c(19)~"LYMPHO",
                            seurat_clusters %in% c(7,8)~"IMM",
                            seurat_clusters %in% c(1,5)~"AST",
                            seurat_clusters %in% c(2)~"OPC",
                            seurat_clusters %in% c(6,20,17,13,16,21)~"INH NEU",
                            seurat_clusters %in% c(11,3,9,12,15)~"EXC NEU")) %>%
  pull(cellid)

data.combined$disease <- data.combined@meta.data %>%
  mutate(disease = case_when(pathology %in% c("control")~"CTRL",
                             T~"MS")) %>%
  pull(disease)

DimPlot(data.combined,reduction = "curated",label = T,group.by = "cellid")
DimPlot(data.combined,reduction = "curated",label = T,group.by = "cellid",split.by = "disease")

# check the ages for the dataset
data.combined@meta.data %>%
  group_by(sample,pathology,disease,sex,age) %>%
  summarise() %>%
  filter(pathology == "control") %>%
  arrange(age)

# load the siganture file
list_sig <- readRDS("../../data/signatures/senescence_pathways.rds")

# wrangling ---------------------------------------------------------------
# add the new classification to the metadata
meta <- data.combined@meta.data %>%
  rownames_to_column("barcodes")

# meta_full <- left_join(meta,LUT,by=c("official_id"))
meta_full <- meta

# add to the original dataset
# data.combined$pathology <- meta_full$pathology

# score the siganture in the UMAP -----------------------------------------
# run the enrichment for the signature. do it on the UMAP using the score siganatures
DefaultAssay(data.combined) <- "RNA"
# x <- "senmayo"

# run the snippet over the whole signatures
lapply(names(list_sig),function(x){
  # extract the dataframe of genes in the signature
  signature.genes.df <- list_sig[[x]]
  
  # pull the genes
  signature.genes <- signature.genes.df %>%
    pull(Genes) %>%
    unique()
  
  # score the module
  data.combined <- AddModuleScore(data.combined,
                                  features = list(signature.genes),
                                  name="signature_score")
  
  # confirm the addition of the score for the module
  # data.combined@meta.data
  df_meta <- data.combined@meta.data %>%
    rownames_to_column("barcode") %>% 
    mutate(signature = x) 
  # mutate(pathology = factor(pathology,levels = c("control cortex","myelinated cortex","demyelinated cortex")))
  
  # save the table with the scores
  df_meta %>% 
    write_tsv(paste0("../../out/table/revision/modules_SENESCENCE_shirmer/120_Module_score_",x,".tsv"))
  
  # save the UMAP coordinates
  df_UMAP <- data.combined@reductions$curated@cell.embeddings %>%
    data.frame() %>%
    rownames_to_column("barcode")
  
  # data2 <- left_join(df_UMAP,df_meta,"barcode")
  # data2_avg <- data2 %>% group_by(cellid) %>% dplyr::select(curated_1, curated_2) %>% summarize_all(mean)
  data2 <- left_join(df_UMAP,df_meta,"barcode")
  data2_avg <- data2 %>% group_by(cellid) %>% dplyr::select(curated_1, curated_2) %>% summarize_all(mean)
  
  # data2 %>%
  #   arrange(signature_score1) %>%
  #   # mutate(gene = "Ptx3") %>%
  #   ggplot() + geom_point(aes(x = curated_1, y = curated_2,col = signature_score1),alpha = 0.5,size = 0.2) +
  #   geom_text_repel(data = data2_avg,aes(x = curated_1, y = curated_2,label = cellid)) +
  #   facet_grid(~pathology) + theme_bw() + 
  #   theme(strip.background = element_blank()) +
  #   # scale_color_gradientn(colours = c("blue","gray","red"))
  #   scale_color_gradientn(colours = viridis::turbo(10))
  # # scale_color_gradientn(colours = c("blue","gray","red"), 
  # #                       values = rescale(c(-0.1,0,0.58)),
  # #                       guide = "colorbar", limits=c(-0.1,0.58))
  # ggsave(paste0("out/image/modules/UMAP_score_",x,".pdf"),width = 12,height = 3)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = curated_1, y = curated_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = curated_1, y = curated_2,label = cellid)) +
    facet_grid(~pathology) + theme_bw() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/120_UMAP_score_",x,".pdf"),width = 10,height = 3)
  ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/120_UMAP_score_",x,".png"),width = 10,height = 3)
  
  # split also by simple ms vs control
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = curated_1, y = curated_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = curated_1, y = curated_2,label = cellid)) +
    facet_grid(~disease) + theme_bw() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/120_UMAP_score_",x,"_2.pdf"),width = 7,height = 3)
  ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/120_UMAP_score_",x,"_2.png"),width = 7,height = 3)
  
  # data2 %>%
  #   arrange(signature_score1) %>%
  #   # mutate(gene = "Ptx3") %>%
  #   ggplot() + geom_point(aes(x = curated_1, y = curated_2,col = signature_score1),alpha = 0.5,size = 0.2) +
  #   geom_text_repel(data = data2_avg,aes(x = curated_1, y = curated_2,label = cellid)) +
  #   facet_wrap(pathology~patient) + theme_bw() + 
  #   theme(strip.background = element_blank()) +
  #   # scale_color_gradientn(colours = c("blue","gray","red"))
  #   scale_color_gradientn(colours = viridis::turbo(10))
  # ggsave(paste0("out/image/modules/UMAP_score_split_",x,".pdf"),width = 12,height = 9)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = curated_1, y = curated_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = curated_1, y = curated_2,label = cellid)) +
    facet_wrap(pathology~sample) + theme_bw() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/120_UMAP_score_split_",x,".pdf"),width = 20,height = 20)
  ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/120_UMAP_score_split_",x,".png"),width = 20,height = 20)
  
  # # plot the score as distribution
  # data2 %>%
  #   # mutate(BraakStage=as.factor(BraakStage)) %>% 
  #   ggplot(aes(x=signature_score1,fill=disease))+geom_density(alpha=0.5)+
  #   facet_grid(origin~cellid,scales = "free")+
  #   theme_bw()+
  #   theme(strip.background = element_blank(), 
  #         panel.border = element_rect(colour = "black", fill = NA))+
  #   scale_fill_manual(values = c("blue","red"))
  # ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE_shirmer/HVG4000/dist_score_",x,"_2.pdf"),width = 12,height =9)
  
  # plot the score as distribution but as ridges
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    ggplot(aes(x=signature_score1,y=disease,fill=disease))+
    ggridges::geom_density_ridges(alpha=0.5)+
    facet_wrap(~cellid,scales = "free")+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    scale_fill_manual(values = c("cyan","red"))
  ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/120_dist_score_ridges_",x,"_2.pdf"),width = 12,height =9)
  
  # plot the score as distribution but as ridges
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    ggplot(aes(x=signature_score1,y=pathology,fill=pathology))+
    ggridges::geom_density_ridges(alpha=0.5)+
    facet_wrap(~cellid,scales = "free")+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))
  # scale_fill_manual(values = c("green","yellow","red"))
  ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/120_dist_score_ridges_",x,".pdf"),width = 12,height =9)
  
  # # plot the distribution of the score per cluster
  # data2 %>%
  #   # mutate(BraakStage=as.factor(BraakStage)) %>% 
  #   arrange(signature_score1) %>%
  #   # mutate(gene = "Ptx3") %>%
  #   mutate(name_sample = paste0(pathology,"_",orig.ident)) %>% 
  #   mutate(rank = rank(pathology,ties.method = "min")) %>% 
  #   mutate(name_sample = fct_reorder(name_sample, rank,.desc = F)) %>% 
  #   
  #   ggplot(aes(x=name_sample,y = signature_score1,fill=pathology)) + 
  #   geom_violin() +
  #   geom_boxplot(width=0.1,outlier.shape = NA) +
  #   facet_wrap(~cellid) + theme_bw() + 
  #   # scale_fill_manual(values = c("green","yellow","red"))+
  #   theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))
  # ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE_shirmer/HVG4000/violin_score_split_",x,".pdf"),width = 25,height = 9)
  
  # plot the distribution of the score per cluster
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot(aes(x=pathology,y = signature_score1,fill=pathology)) + 
    geom_violin() +
    geom_boxplot(width=0.1,outlier.shape = NA) +
    facet_wrap(~cellid) + theme_bw() + 
    # scale_fill_manual(values = c("green","yellow","red"))+
    theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))
  ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/120_violin_score_",x,".pdf"),width = 9,height = 9)
})

# tailored plotting -------------------------------------------------------
x <- "senmayo"

signature.genes.df <- list_sig[[x]]

# pull the genes
signature.genes <- signature.genes.df %>%
  pull(Genes) %>%
  unique()

# score the module
data.combined <- AddModuleScore(data.combined,
                                features = list(signature.genes),
                                name="signature_score")

# confirm the addition of the score for the module
# data.combined@meta.data
df_meta <- data.combined@meta.data %>%
  rownames_to_column("barcode") %>% 
  mutate(signature = x) 
# mutate(pathology = factor(pathology,levels = c("control cortex","myelinated cortex","demyelinated cortex")))

# save the UMAP coordinates
df_UMAP <- data.combined@reductions$curated@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")

# data2 <- left_join(df_UMAP,df_meta,"barcode")
# data2_avg <- data2 %>% group_by(cellid) %>% dplyr::select(curated_1, curated_2) %>% summarize_all(mean)
data2 <- left_join(df_UMAP,df_meta,"barcode")
data2_avg <- data2 %>% group_by(cellid) %>% dplyr::select(curated_1, curated_2) %>% summarize_all(mean)

# data2 %>%
#   arrange(signature_score1) %>%
#   # mutate(gene = "Ptx3") %>%
#   ggplot() + geom_point(aes(x = curated_1, y = curated_2,col = signature_score1),alpha = 0.5,size = 0.2) +
#   geom_text_repel(data = data2_avg,aes(x = curated_1, y = curated_2,label = cellid)) +
#   facet_grid(~pathology) + theme_bw() + 
#   theme(strip.background = element_blank()) +
#   # scale_color_gradientn(colours = c("blue","gray","red"))
#   scale_color_gradientn(colours = viridis::turbo(10))
# # scale_color_gradientn(colours = c("blue","gray","red"), 
# #                       values = rescale(c(-0.1,0,0.58)),
# #                       guide = "colorbar", limits=c(-0.1,0.58))
# ggsave(paste0("out/image/modules/UMAP_score_",x,".pdf"),width = 12,height = 3)

data2 %>%
  arrange(signature_score1) %>%
  mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active"))) %>% 
  # mutate(gene = "Ptx3") %>%
  ggplot() + geom_point(aes(x = curated_1, y = curated_2,col = signature_score1),alpha = 0.5,size = 0.2) +
  # geom_text_repel(data = data2_avg,aes(x = curated_1, y = curated_2,label = cellid)) +
  facet_wrap(~pathology,nrow = 3) + theme_void() + 
  theme(strip.background = element_blank()) +
  scale_color_gradientn("sig score",
                        colours = c("gray","orange","red"),
                        oob = scales::squish,limits = c(0,0.25))
# scale_color_gradientn(colours = c("gray","yellow","red"))
# scale_color_gradientn(colours = viridis::turbo(10))
ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/119_UMAP_score_",x,"_tailored1.pdf"),width = 4,height = 9)
ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/119_UMAP_score_",x,"_tailored1.png"),width = 4,height = 9,bg="white")

# data2 %>%
#   arrange(signature_score1) %>%
#   mutate(pathology = factor(pathology,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CA","WM_CI","WM_Core"))) %>% 
#   # mutate(gene = "Ptx3") %>%
#   ggplot() + geom_point(aes(x = curated_1, y = curated_2,col = signature_score1),alpha = 0.5,size = 0.2) +
#   geom_text_repel(data = data2_avg,aes(x = curated_1, y = curated_2,label = cellid)) +
#   facet_wrap(~pathology,nrow = 3) + theme_void() + 
#   theme(strip.background = element_blank()) +
#   scale_color_gradientn("sig score",
#                         colours = viridis::turbo(10),
#                         values = rescale(c(-0.1,0,0.3)),
#                         oob = scales::squish,limits = c(-0.1,0.3))

data2 %>%
  arrange(signature_score1) %>%
  mutate(pathology = factor(pathology,levels = c("control","chronic_inactive","chronic_active"))) %>% 
  # mutate(gene = "Ptx3") %>%
  ggplot() + geom_point(aes(x = curated_1, y = curated_2,col = signature_score1),alpha = 0.5,size = 0.2) +
  # geom_text_repel(data = data2_avg,aes(x = curated_1, y = curated_2,label = cellid)) +
  facet_wrap(~pathology,nrow = 3) + theme_void() + 
  theme(strip.background = element_blank()) +
  scale_color_gradientn("sig score",
                        colours = viridis::turbo(10),limits = c(-0.1,0.3),oob = scales::squish)
ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/119_UMAP_score_",x,"_tailored2.pdf"),width = 4,height = 9)
ggsave(paste0("../../out/image/revision/modules_SENESCENCE_shirmer/119_UMAP_score_",x,"_tailored2.png"),width = 4,height = 9,bg="white")

# run Upset plot for the signatures ---------------------------------------
# load the siganture file
# list_sig <- readRDS("../../data/signatures/senescence_pathways.rds")

# shortlist the signatures and rename them
# list_sig_shortlist <- list_sig %>%
#   bind_rows(.id = "signature") %>%
#   filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>%
#   mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
#                                 signature == "Inhibits"~"CellAge Inhibits",
#                                 signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
#                                 signature == "senmayo"~"Senmayo")) %>%
#   split(f = .$signature2) %>%
#   lapply(function(x){
#     x %>% pull(Genes) %>%
#       unique()
#   })


# library(UpSetR)
# pdf("../../out/image/revision/modules_SENESCENCE_shirmer/120_upset_signatures_shortlist.pdf",width = 5,height = 4,onefile = F)
# upset(fromList(list_sig_shortlist), order.by = "freq",nsets = 10,nintersects = 100)
# dev.off()

# pull the intersections
# str(list_sig_shortlist)
# 
# df1 <- lapply(list_sig_shortlist,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# head(df1)
# df2 <- data.frame(gene=unique(unlist(list_sig_shortlist)))
# 
# head(df2)

# now loop through each individual gene and pick the list of all the intersections they belong to
# df_int <- lapply(df2$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1 %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
#   
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()

# head(df_int,n=20)

# save the intersection table
# df_int %>%
#   write_tsv("../../out/table/revision/modules_SENESCENCE_shirmer/120_upset_signatures_shortlist_intersection.tsv")

# confirm the data and the list are congruent
# df_int %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))

# it is in keeping with the expected result
