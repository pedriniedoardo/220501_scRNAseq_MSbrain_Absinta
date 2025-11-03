# AIM ---------------------------------------------------------------------
# test bulk deconvolution using BayesPrism
# to make the computation faster, use the GEP approach rather than the count.matrix one
# in this case use the VAS subcluster profile

# libraries ---------------------------------------------------------------
library(Seurat)
library(BayesPrism)
library(tidyverse)

# read in the data --------------------------------------------------------
# # read in the data from the sc dataset
# scobj <- readRDS("../../out/object/130_VAS_subcluster_HarmonySample.rds")
# DimPlot(scobj,raster = T,group.by = "RNA_snn_res.0.4",label = T)
# 
# # remember that according to the signature scoring we have proposed the the main endothelial population is composed of cells from cluster:
# # 0: capillaries
# # 5: veins
# # 9: arterial
# # all the other clusters are likely contaminats
# 
# # simple aggregation of the data by cell type
# cts_sample_CellidSmall <- AggregateExpression(object = scobj,
#                                               group.by = c("RNA_snn_res.0.4"),
#                                               assays = 'RNA',
#                                               slot = "counts",
#                                               return.seurat = FALSE)
# 
# # extract the table of counts
# ref_GEP <- cts_sample_CellidSmall$RNA
# 
# # save the aggregated matrix
# saveRDS(ref_GEP,"../../out/object/207_ref_GEP_VAS_subcluster.rds")

# read in the reference bulk dataset
ref_GEP <- readRDS("../../out/object/207_ref_GEP_VAS_subcluster.rds")

# read in the sample bulk test data. 
# we need at least two samples
# dds <- readRDS(file = "../../out/object/201_dds_all.rds")
dds <- readRDS("../../out/object/201_dds_all_update.rds")
bulk_counts <- counts(dds, normalized = F)

# Check the dimensions of our simulated data
dim(ref_GEP)
dim(bulk_counts)

#  Step 2: Prepare Data for BayesPrism ------------------------------------
# BayesPrism requires the data in specific formats. We need to extract the raw counts and cell labels from our Seurat object and ensure the gene names match between our reference and bulk data.

# Get cell type labels. These will be our 'cell.type.labels'.
cell_types <- colnames(ref_GEP) %>% unname()

# Get cell state labels. For a simple case, these can be the same as cell types.
# In more complex analyses, these could be cluster IDs or subtypes.
# https://github.com/Danko-Lab/BayesPrism/issues/66
cell_states <- colnames(ref_GEP) %>% unname()

# Filter and Align Genes
# reshape the data to be accepted by the tool
# generally data are reporeted as feature X sample, but in the tool they are handled as sample X feature.
# we can input sparse matrices: https://github.com/Danko-Lab/BayesPrism/issues/58
GEP.dat <- t(ref_GEP)
GEP.dat[1:5,1:5]
dim(GEP.dat)

bk.dat <- t(bulk_counts)
bk.dat[,1:5]
dim(bk.dat)

# Run the reccommended filtering on the genes in this step for a real data analysis
# Filter outlier genes from ref data
# Next, we remove the genes from selected groups. Note that when sex is not identical between the reference and mixture, we recommend excluding genes from chrX and chrY. We also remove lowly transcribed genes, as the measurement of transcription of these genes tend to be noise-prone. Removal of these genes can also speed up computation.

GEP.dat.filtered <- cleanup.genes(input = GEP.dat,
                                  input.type = "GEP",
                                  species = "hs",
                                  gene.group = c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"))

dim(GEP.dat.filtered)
dim(GEP.dat)

# Next, we check the concordance of gene expression for different types of genes. As bulk and single cell data are usually collected by different experimental protocols, they may have different sensitivity to different types of genes.
# note this function only works for human data. For other species, you are advised to make plots by yourself.
plot.bulk.vs.sc(sc.input = GEP.dat.filtered,
                bulk.input = bk.dat
                #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)

# Subset protein coding genes.
GEP.dat.filtered.pc <-  select.gene.type(GEP.dat.filtered,
                                         gene.type = "protein_coding")
str(GEP.dat.filtered.pc)
dim(GEP.dat.filtered.pc)
dim(GEP.dat.filtered)
dim(GEP.dat)

# potentially further filter the dataset to include only known markers genes per cell type

# Step 3: Run BayesPrism Deconvolution ------------------------------------
# With the data properly formatted and aligned, we can now run the deconvolution.

# Construct the BayesPrism Object
# We combine the reference and cell labels into a 'prism' object.
my_prism <- new.prism(
  reference = GEP.dat.filtered.pc,
  mixture = bk.dat,
  key = NULL, # No specific tumor/normal key needed for this example
  cell.type.labels = cell_types,
  cell.state.labels = cell_states,
  input.type = "GEP"
)

# Run the Deconvolution Algorithm
# This is the core computational step.
# We'll use 8 cores for this example. Adjust 'n.cores' as needed.
bp_results <- run.prism(prism = my_prism, n.cores = 8)

# You can look at the returned object to see what it contains
bp_results

# saveRDS(bp_results,"../out/object/bp_results_rawCount_pBulk_treat_100_ifnb.rds")
saveRDS(bp_results,"../../out/object/207_bp_results_rawCount_VAS_sublsuter_GEP.rds")

# Step 4: Extract and Calculate Deconvoluted Reads ------------------------
# this can be extracted from the object as follows
# x <- "Endothelial"

list_bulk <- lapply(unique(cell_types), function(x){
  mat <- get.exp(bp=bp_results,
                 state.or.type="type",
                 cell.name=x)
  
  data <- mat %>%
    t()
  
  return(data)
}) %>%
  setNames(unique(cell_types))

# generate integer estimates
list_bulk_integer <- lapply(unique(cell_types), function(x){
  mat <- get.exp(bp=bp_results,
                 state.or.type="type",
                 cell.name=x)
  
  # generate integer level estimates
  # https://github.com/Danko-Lab/BayesPrism/issues/121
  # round up the Z matrix to the nearest integer and use it as input for DESeq2 (without any normalization step).
  data <- mat %>%
    t() %>%
    round()
  return(data)
}) %>%
  setNames(unique(cell_types))

# notice that the number of genes matches the one present in the reference
lapply(list_bulk, function(x){
  dim(x)
})

dim(GEP.dat.filtered.pc)

# get the fraction estimated
theta <- get.fraction(bp=bp_results,
                      which.theta="final",
                      state.or.type="type")
theta

# save the prop 
df_theta <- theta %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(names_to = "expertAnno.l1",values_to = "prop",-sample)

# save the table
write_tsv(df_theta,"../../out/table/207_theta_VAS_subcluster_GEP.tsv")

# plot the porportion calculated based on the deconvolution
df_theta %>%
  ggplot(aes(x=expertAnno.l1,y=prop))+geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),alpha=0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave(filename = "../../out/image/207_boxplot_theta_VAS_subcluster_GEP.pdf",width = 5,height = 5)  
