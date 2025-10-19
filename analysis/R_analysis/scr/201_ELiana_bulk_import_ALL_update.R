# AIM ---------------------------------------------------------------------
# import the raw data shared by Eliana, relative to her bulk dataset
# this is the analysis for the MS vs CTRL samples
# this is the updata version of the analysis

# libraries ---------------------------------------------------------------
library("tidyverse")
library("DESeq2")

# read in the expression data ---------------------------------------------
# read in the raw table of counts
raw_counts <- readRDS("../../data/eliana_bulkRNAseq/eliana_mat.rds")
head(raw_counts)

# metadata ----------------------------------------------------------------
# load the updated metadata. according to Eliana they have removed sample MS_96 and kept sample CTR_33 but they switched the gender
LUT_samples <- read_csv("../../data/eliana_bulkRNAseq/LUT_samples_update.csv") %>%
  column_to_rownames()

# check the balance of the covariates
# in the first comparison
LUT_samples %>%
  group_by(Group,sex) %>%
  summarise(n = n())

# in the second comparison
LUT_samples %>%
  group_by(Treatment,sex) %>%
  summarise(n = n())

# match the meta with the count table
mat_exp <- raw_counts %>% 
  .[,LUT_samples$`HUGE-ID`]

# define the model --------------------------------------------------------
# clone <- coldata$clone
gender <- LUT_samples$sex
# notice that for this analysis the treatment variable cannot be included as it is colinar with the disease variable
# treat <- factor(LUT_samples$Treatment)
disease <- factor(LUT_samples$Group,levels = c("CTRL","MS"))

# build the design
design <- model.matrix(~ gender + disease)
colnames(design)[1] <- c("intercept")

saveRDS(design,file = "../../out/object/201_design_all_update.rds")

# build the object --------------------------------------------------------
# is keeping only the objext in the lut_sample
dds <- DESeqDataSetFromMatrix(countData = mat_exp,
                              colData = LUT_samples,
                              design = design)

saveRDS(dds,file = "../../out/object/201_dds_all_update.rds")
