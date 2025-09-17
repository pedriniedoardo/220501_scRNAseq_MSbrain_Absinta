# libraries ---------------------------------------------------------------
library(tidyverse)
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)

# read in the data --------------------------------------------------------
# scobj <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegration_AllSoupX_test.rds")
# 
# # add the cell_id based on martina's annotation
# scobj$cell_id <- scobj@meta.data %>%
#   mutate(cell_id = case_when(seurat_clusters %in% c(0,4,14)~"OLIGO",
#                              seurat_clusters %in% c(9)~"OPC",
#                              seurat_clusters %in% c(3,12)~"ASTRO",
#                              seurat_clusters %in% c(5)~"IMMUNE",
#                              seurat_clusters %in% c(13)~"LYM",
#                              seurat_clusters %in% c(11)~"VAS",
#                              seurat_clusters %in% c(1, 2, 10,6)~"EXC NEU",
#                              seurat_clusters %in% c(7,8)~"INH NEU")) %>% 
#   pull(cell_id)
# 
# 
# meta <- scobj@meta.data %>%
#   rownames_to_column("barcode")
# write_tsv(meta,"../../out/table/data.combined_WM_CX_harmonySkipIntegration_AllSoupX_test_propeller.tsv")
# 
meta_ref <- read_tsv(file = "../../out/table/data.combined_WM_CX_harmonySkipIntegration_AllSoupX_test_propeller.tsv")


# test 1 ------------------------------------------------------------------
# confirm the numbers from the tissue dataset
meta_ref %>% 
  group_by(orig.ident,origin,cell_id) %>% 
  summarise(n = n())

# brain MS vs control -----------------------------------------------------
# run it on the predictred l2, not the robust one to avoid the missing values
# Run propeller testing for cell type proportion differences between the groups
out1 <- propeller(clusters = meta_ref$cell_id,
                 sample = meta_ref$orig.ident,
                 group = meta_ref$origin)

out1 %>%
  rownames_to_column("cluster") %>%
  write_tsv("../../out/table/propeller_out_tissue.tsv")


# plotting
# plotCellTypeProps(clusters = meta_ref$assigned.predicted.id.l1,
#                   sample = meta_ref$orig.ident.cca)
df_summary <- meta_ref %>% 
  group_by(cell_id,
           orig.ident,
           origin) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(orig.ident) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot
df_summary %>%
  ggplot(aes(x=origin,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1))+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/propeller_plot.pdf",width = 6,height = 6)

# df_summary %>%
#   ggplot(aes(x=origin,y=prop))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1))+
#   geom_text_repel(aes(label = orig.ident,x=origin,y=prop))+
#   facet_wrap(~cell_id,scales = "free")+
#   theme_bw()+
#   theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# # ggsave("out/image/propeller_plot2.pdf",width = 20,height = 15)

df_summary %>%
  # mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>% 
  ungroup() %>% 
  group_by(cell_id,origin) %>% 
  summarise(n = sum(n)) %>% 
  ungroup() %>% 
  group_by(origin) %>% 
  mutate(tot = sum(n),
         prop = n/tot) %>% 
  ggplot(aes(x=origin,y=prop,fill=cell_id))+
  geom_bar(stat = "identity", position="fill") +
  ggtitle("SCType predictions") +
  theme_cowplot() +
  labs(x = "SampleID", y = "fraction of cells") 
ggsave("../../out/image/propeller_plot_stacked.pdf",width = 4,height = 4)
