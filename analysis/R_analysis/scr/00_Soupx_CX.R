# libraries ---------------------------------------------------------------
library(tidyverse)
library(SoupX)

# READ IN DATA ------------------------------------------------------------
id_sample <- dir("../../data/cellranger_minimal/sample_cellranger6.1/")

# load the LUT
# LUT <- read_csv("data/clinical_data.csv") %>% 
#   mutate(sample_id = str_pad(ID,width = 2,side = "left",pad = "0")) %>% 
#   mutate(sample_id = paste0(sample_id,"_cr_61"))

# do the preprocessing over all the dataset and save the output matrices
# x <- "MA8483-8_all"
# file <- paste0("../../data/cellranger_minimal/sample_cellranger6.1/",x,"/outs/")
# test <- load10X(dataDir = file)
# test2 <- autoEstCont(test)
# test3 <- adjustCounts(test2,roundToInt = T)

lapply(id_sample,function(x){
  # define the location of the output of cellranger
  file <- paste0("../../data/cellranger_minimal/sample_cellranger6.1/",x,"/outs/")
  # Load data and estimate soup profile
  sc <- load10X(file)
  # Estimate rho
  sc <- autoEstCont(sc)
  # Clean the data
  out <- adjustCounts(sc,roundToInt = T)
  # save the data
  DropletUtils:::write10xCounts(paste0("../../data/SoupX_default/",x), out)
})
