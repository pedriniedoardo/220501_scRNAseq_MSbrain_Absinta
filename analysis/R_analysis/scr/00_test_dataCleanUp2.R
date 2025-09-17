# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggraph)
library(igraph)
library(data.tree)
library(HGNChelper)

# read in the data --------------------------------------------------------
scobj_all <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_harmonySkipIntegration_AllSoupX_test.rds")

# test on a small subset
scobj_subset <- subset(scobj_all,orig.ident == "s10")

# confirm the identity of the dataset
DimPlot(scobj_all,raster=T,label = T)

# processing --------------------------------------------------------------
#clean dataset
# all20 <- scobj_subset
all20 <- scobj_all
pv.cutoff.2.use <- 0.1336 # edit here
all20.var.genes <- VariableFeatures(all20)
all20.scale.data <- GetAssayData(all20,slot="scale.data")
all20.scale.data <- all20.scale.data[all20.var.genes,]
all20.clusters <- all20$seurat_clusters
all20.clusters.unique <- unique(sort(all20.clusters))
to.keep <- NULL

for(i in 1:length(all20.clusters.unique)) {
  print(i)
  temp1 <- all20.scale.data[,all20.clusters==all20.clusters.unique[i]]
  temp2 <- cor(temp1)
  dim(temp2)
  temp3 <- NULL
  for(j in 1:dim(temp2)[2]) {
    temp4 <- temp2[,j][-c(j)]
    temp5 <- mean(temp4)
    temp3 <- c(temp3,temp5)
  }
  temp6 <- NULL
  temp100 <- NULL
  for(k in 1:length(temp3)) {
    temp7 <- temp3[k]
    temp8 <- mean(temp3[-c(k)])
    temp9 <- sd(temp3[-c(k)])
    temp10 <- (temp7-temp8)/temp9
    temp11 <- 1-pnorm(abs(temp10))
    temp6 <- c(temp6,temp11)
    temp100 <- c(temp100,temp10)
  }
  names(temp6) <- dimnames(temp1)[[2]]
  names(temp100) <- dimnames(temp1)[[2]]
  temp12 <- names(temp6[temp6<pv.cutoff.2.use])
  temp12 <- temp100[temp12]
  temp12 <- temp12[temp12<0]
  temp12 <- setdiff(dimnames(temp1)[[2]],names(temp12))
  to.keep <- c(to.keep,temp12)
}
temp13 <- rep(1,length(to.keep))
names(temp13) <- to.keep
temp14 <- rep(0,(dim(all20.scale.data)[2]-length(temp13)))
names(temp14) <- setdiff(dimnames(all20.scale.data)[[2]],to.keep)
temp15 <- c(temp13,temp14)
temp16 <- temp15[dimnames(all20.scale.data)[[2]]]
all20@meta.data$curation.pass <- temp16
all20
z.score.curated.by.cluster.all20 <- subset(all20, subset = curation.pass > 0)
z.score.curated.by.cluster.all20
z.score.curated.by.cluster.all20 <- RunUMAP(z.score.curated.by.cluster.all20, dims = 1:10,reduction.name="curated")
DimPlot(z.score.curated.by.cluster.all20, reduction = "curated")

# save the final object ---------------------------------------------------
saveRDS(z.score.curated.by.cluster.all20,"../../out/object/00_test_dataCleanUP_test.rds")

