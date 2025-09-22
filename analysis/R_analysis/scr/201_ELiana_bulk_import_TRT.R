# AIM ---------------------------------------------------------------------
# import the raw data shared by Eliana, relative to her bulk dataset
# this is the analysis for the Naive vs Treated samples, therefore I will need to trimm the samples

# libraries ---------------------------------------------------------------
library("tidyverse")
library("DESeq2")

# read in the expression data ---------------------------------------------
# read in the raw table of counts
raw_counts <- readRDS("../../data/eliana_bulkRNAseq/eliana_mat.rds")
head(raw_counts)

# metadata ----------------------------------------------------------------
# build the annotation besed on the sample metadata
LUT_samples <- readRDS("../../data/eliana_bulkRNAseq/eliana_metadata.rds") %>%
  # filter only the samples in Naive and TRT
  dplyr::filter(Treatment %in% c("TRT","Naive"))

# check the balance of the covariates
LUT_samples %>%
  group_by(Treatment,sex) %>%
  summarise(n = n())

# match the meta with the count table
mat_exp <- raw_counts %>% 
  .[,LUT_samples$`HUGE-ID`]

# define the model --------------------------------------------------------
# clone <- coldata$clone
gender <- LUT_samples$sex
treat <- factor(LUT_samples$Treatment)
# disease <- factor(LUT_samples$Group,levels = c("CTRL","MS"))

# build the design
design <- model.matrix(~ gender + treat)
colnames(design)[1] <- c("intercept")

saveRDS(design,file = "../../out/object/201_design_trt.rds")

# build the object --------------------------------------------------------
# is keeping only the objext in the lut_sample
dds <- DESeqDataSetFromMatrix(countData = mat_exp,
                              colData = LUT_samples,
                              design = design)

saveRDS(dds,file = "../../out/object/201_dds_trt.rds")
