# This vignette shows how to apply CellChat to identify major signaling changes as well as conserved and context-specific signaling by joint manifCSF learning and quantitative contrasts of multiple cell-cell communication networks. We showcase CellChat’s diverse functionalities by applying it to a scRNA-seq data on cells from two biological conditions: nonlesional (NL, normal) and lesional (LS, diseased) human skin from patients with atopic dermatitis. These two datasets (conditions) have the same cell population compositions after joint clustering. If there are slightly or vastly different cell population compositions between different datasets, please check out another related tutorial (Comparison analysis of multiple datasets with different cell type compositions).

# CellChat employs a top-down approach, i.e., starting with the big picture and then refining it in a greater detail on the signaling mechanisms, to identify signaling changes at different levels, including both general principles of cell-cell communication and dysfunctional cell populations/signaling pathways/ligand-receptors.

# Load the required libraries ---------------------------------------------
library(CellChat)
library(patchwork)
library(circlize)

# Create a directory to save figures --------------------------------------
# data.dir <- './cellchat_comparison_out'
# dir.create(data.dir)
# setwd(data.dir)

# Load CellChat object of each dataset and then merge together ------------
# USERS need to run CellChat on each dataset seperately and then merge different CellChat objects together. Please do updateCellChat if you have CellChat objects that are obtained using the earlier version (< 0.5.0).

# cellchat.NL <- readRDS(url("https://ndownloader.figshare.com/files/25954199"))
# saveRDS(cellchat.NL,file = "../data/cellchat_humanSkin_NL.rds")
# cellchat.LS <- readRDS(url("https://ndownloader.figshare.com/files/25956518"))
# saveRDS(cellchat.LS,file = "../data/cellchat_humanSkin_LS.rds")
cellchat.CTRL <- readRDS("../../out/object/revision/126_cellchat_GroupCTRLCellID_full.rds")
cellchat.SEN <- readRDS("../../out/object/revision/126_cellchat_GroupSENCellID_full.rds")

# explore the entity of the rds objects
cellchat.CTRL
cellchat.SEN

# these are cellchat objects that can be generated usign:
# cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
# for more details see the tutorial at C:\Users\pedri\OneDrive\Documents\training\bioinformatic\R_package_CellChat\CellChat\Inference_and_analysis_of_cell-cell_communication_using_CellChat\sample_workflow.R

# put the objects in a list
object.list <- list(CTRL = cellchat.CTRL, SEN = cellchat.SEN)

# merge the objects
# notice the same barcode cannot be used in the same datasets
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

# print the content of the object
cellchat

# Part I: Predict general principles of cell-cell communication -----------
# CellChat starts with the big picture to predict general principles of cell-cell communication. When comparing cell-cell communication among multiple biological conditions, it can answer the following biological questions:
# Whether the cell-cell communication is enhanced or not
# The interaction between which cell types is significantly changed
# How the major sources and targets change from one condition to another

# Compare the total number of interactions and interaction strength -------
# To answer on question on whether the cell-cell communication is enhanced or not, CellChat compares the the total number of interactions and interaction strength of the inferred cell-cell communication networks from different biological conditions.

# measure: "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

gg1 + gg2
ggsave("../../out/image/revision/128_comparison_01_barplot_interaction_cellID.pdf",width = 10,height = 5)

# group: a vector giving the groups of different datasets to define colors of the bar plot. Default: only one group and a single color
# how is the group argument affectiong the output?
# compareInteractions(cellchat, show.legend = F)
# compareInteractions(cellchat, show.legend = F, group = c(1,2))

# Compare the number of interactions and interaction strength amon --------
# Compare the number of interactions and interaction strength among different cell populations
# To identify the interaction between which cell populations showing significant changes, CellChat compares the number of interactions and interaction strength among different cell populations.

# Differential number of interactions or interaction strength amon --------
# Differential number of interactions or interaction strength among different cell populations
# The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.

# pdf("02_differential_interaction.pdf",width = 10,height = 10)
# par(mfrow = c(1,2), xpd=TRUE)
# lapply(c("count","weight"),function(i){
#   netVisual_diffInteraction(cellchat, weight.scale = T, measure = i)
# })
# dev.off()

# pdf("../../out/image/revision/128_comparison_02_differential_interaction_cellID.pdf",width = 20,height = 10)
# par(mfrow = c(1,2), xpd=TRUE)
# list(netVisual_diffInteraction_fix(cellchat, weight.scale = T),
#      netVisual_diffInteraction_fix(cellchat, weight.scale = T, measure = "weight"))
# dev.off()

pdf("../../out/image/revision/128_comparison_02_differential_interaction_cellID.pdf",width = 20,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
list(netVisual_diffInteraction(cellchat, weight.scale = T),
     netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight"))
dev.off()

# try to reduce the scale
pdf("../../out/image/revision/128_comparison_02_differential_interaction_cellID_high.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
list(netVisual_diffInteraction(cellchat, weight.scale = T),
     netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight"))
dev.off()

netVisual_diffInteraction(cellchat, weight.scale = T,edge.width.max = 20,label.edge = T)

pdf("../../out/image/revision/128_comparison_02_differential_interaction_cellID_high2.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
list(netVisual_diffInteraction(cellchat, weight.scale = T,edge.width.max = 20,vertex.label.cex = 0.01))
dev.off()

# netVisual_diffInteraction(cellchat, weight.scale = T)
# netVisual_diffInteraction(cellchat, weight.scale = T,sources.use = "Endo",targets.use = "Endo")

pdf("../../out/image/revision/128_comparison_02_differential_interaction_UP_cellID.pdf",width = 20,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
list(fun_test(cellchat, weight.scale = T,filter_edges = "up", measure = "count"),
     fun_test(cellchat, weight.scale = T,filter_edges = "up", measure = "weight"))
dev.off()

pdf("../../out/image/revision/128_comparison_02_differential_interaction_DOWN_cellID.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
list(fun_test(cellchat, weight.scale = T,filter_edges = "down"),
     fun_test(cellchat, weight.scale = T,filter_edges = "down", measure = "weight"))
dev.off()

# cellchat@net$AC$prob %>%
#   as.data.frame() %>%
#   dplyr::select(contains("NCAM1_NCAM2"))
# cellchat@netP$AC$prob %>%
#   as.data.frame() %>%
#   dplyr::select(contains("NCAM"))

#
# cellchat@net$CCA$prob %>%
#   as.data.frame() %>%
#   dplyr::select(contains("NCAM1_NCAM2"))
# cellchat@netP$CCA$prob %>%
#   as.data.frame() %>%
#   dplyr::select(contains("NCAM"))

# cellchat.SEN@net$prob %>%
#   as.data.frame() %>%
#   dplyr::select(contains("NCAM1_NCAM2"))
#
# cellchat.CTRL@net$prob %>%
#   as.data.frame() %>%
#   dplyr::select(contains("NCAM1_NCAM2"))

# computeAveExpr(cellchat.CTRL, features = c("NCAM1","NCAM2"), type =  "triMean")
# computeAveExpr(cellchat.SEN, features = c("NCAM1","NCAM2"), type =  "triMean")

# We can also show differential number of interactions or interaction strength in a greater details using a heatmap. The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of row of values (outgoing signaling). In the colorbar, red (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.

gg1 <- netVisual_heatmap(cellchat,height = 5,comparison = c(1,2))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",height = 5)

pdf("../../out/image/revision/128_comparison_02_differential_interaction_heatmap_cellID.pdf",width = 10,height = 5)
gg1 + gg2
dev.off()

# To better control the node size and edge weights of the inferred networks across different datasets, we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
weight.max

pdf("../../out/image/revision/128_comparison_03_interaction_individual_number_cellID.pdf",width = 20,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
lapply(1:length(object.list),function(i){
  netVisual_circle(object.list[[i]]@net$count,
                   weight.scale = T,
                   label.edge= F,
                   edge.weight.max = weight.max[2],
                   edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
})
dev.off()

# same as above but for strength
pdf("../../out/image/revision/128_comparison_03_interaction_individual_strenght_cellID.pdf",width = 20,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
lapply(1:length(object.list),function(i){
  netVisual_circle(object.list[[i]]@net$weight,
                   weight.scale = T,
                   label.edge= F,
                   # edge.weight.max = weight.max[2],
                   edge.width.max = 12,
                   title.name = paste0("Strenght of interactions - ", names(object.list)[i]))
})
dev.off()

# Compare the major sources and targets in 2D space -----------------------
# Comparing the outgoing and incoming interaction strength in 2D space allows ready identification of the cell populations with significant changes in sending or receiving signals between different datasets.

num.link <- sapply(object.list, function(x){
  rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)
})
num.link

# control the dot size in the different datasets
weight.MinMax <- c(min(num.link), max(num.link))
weight.MinMax

gg <- lapply(1:length(object.list),function(i){
  netAnalysis_signalingRole_scatter(object.list[[i]],
                                    title = names(object.list)[i],
                                    weight.MinMax = weight.MinMax)
})

pdf("../../out/image/revision/128_comparison_04_compare_source_target_cellID.pdf",width = 10,height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

# From the scatter plot, we can see that Inflam.DC and cDC1 emerge as one of the major source and targets in LS compared to NL. Fibroblast populations also become the major sources in LS.

# Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and LS. ## Identify signaling changes associated with one cell group

# Visualizing differential outgoing and incoming signaling changes from NL to LS
# gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Inflam. DC", signaling.exclude = "MIF")
# gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))
id <- c("AST","IMM","OLIGO","EXC NEU","INH NEU","VAS","OPC","EPENDYMA","LYM")

list_gg <- lapply(id, function(x){
  netAnalysis_signalingChanges_scatter(cellchat, idents.use = x)
})

# netAnalysis_signalingChanges_scatter(cellchat, idents.use = "clu_14")

pdf("../../out/image/revision/128_comparison_05_compare_signalling_cellID.pdf",width = 30,height = 15)
patchwork::wrap_plots(plots = list_gg)
dev.off()

id2 <- c("IMM")

list_gg2 <- lapply(id2, function(x){
  netAnalysis_signalingChanges_scatter(cellchat, idents.use = x)
})

# netAnalysis_signalingChanges_scatter(cellchat, idents.use = "clu_14")

pdf("../../out/image/revision/128_comparison_05_compare_signalling_cellID2.pdf",width = 12,height = 10)
patchwork::wrap_plots(plots = list_gg2)
dev.off()

# Part II: Identify the conserved and context-specific signaling p --------
# CellChat then can identify signaling networks with larger (or less) difference, signaling groups, and the conserved and context-specific signaling pathways based on their cell-cell communication networks among multiple biological conditions.

# Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
# CellChat performs joint manifCSF learning and classification of the inferred communication networks based on their functional and topological similarity. NB: Such analysis is applicable to more than two datasets.

# Functional similarity: High degree of functional similarity indicates major senders and receivers are similar, and it can be interpreted as the two signaling pathways or two ligand-receptor pairs exhibit similar and/or redundant roles. NB: Functional similarity analysis is not applicable to multiple datsets with different cell type composition.

# Structural similarity: A structural similarity was used to compare their signaling network structure, without considering the similarity of senders and receivers. NB: Structural similarity analysis is applicable to multiple datsets with the same cell type composition or the vastly different cell type composition.

#  Here we can run the manifCSF and classification learning analysis based on the functional similarity because the two datasets have the the same cell type composition.

# Identify signaling groups based on their functional similarity


# issue -------------------------------------------------------------------
# try the fix from the internet fix
# the fixing produce a different plot from the tutorial.
# source("../00_functional.R")
# cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")
#
# # Visualization in 2D-space
# netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)


# # Visualization in 2D-space
# pdf("05_functional_similarity.pdf",width = 10,height = 10)
# netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
# dev.off()

# -------------------------------------------------------------------------

# Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural",do.parallel = F)

# Visualization in 2D-space
pdf("../../out/image/revision/128_comparison_05_structural_similarity_cellID.pdf",width = 10,height = 9)
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()

pdf("../../out/image/revision/128_comparison_05_structural_similarity_grid_cellID.pdf",width = 10,height = 9)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()

# Compute and visualize the pathway distance in the learned joint  --------

# Compute and visualize the pathway distance in the learned joint manifCSF
# We can identify the signaling networks with larger (or less) difference based on their Euclidean distance in the shared two-dimensions space. Larger distance implies larger difference of the communication networks between two datasets in terms of either functional or structure similarity. NB: We only compute the distance of overlapped signaling pathways between two datasets. Those signaling pathways that are only identified in one dataset are not considered here. If there are more than three datasets, one can do pairwise comparisons by defining comparison in the function rankSimilarity.

pdf("../../out/image/revision/128_comparison_06_rank_similarity_cellID.pdf",width = 5,height = 5)
rankSimilarity(cellchat, type = "structural")
dev.off()

# Identify and visualize the conserved and context-specific signal --------

# Identify and visualize the conserved and context-specific signaling pathways
# By comparing the information flow/interaction strengh of each signaling pathway, we can identify signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv) increase, by change their information flow at one condition as compared to another condition.

# Compare the overall information flow of each signaling pathway
# We can identify the conserved and context-specific signaling pathways by simply comparing the information flow for each signaling pathway, which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network).

# This bar graph can be plotted in a stacked mode or not. Significant signaling pathways were ranked based on differences in the overall information flow within the inferred networks between NL and LS skin. The top signaling pathways colored red are enriched in NL skin, and these colored green were enriched in the LS skin.

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

pdf("../../out/image/revision/128_comparison_07_signaling_pathway_comparison_cellID.pdf",width = 10,height = 8)
gg1 + gg2
dev.off()

# levels(cellchat@idents$joint)
# rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,sources.use = 2)
# rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,sources.use = 2)
# rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,targets.use = 1)

# Compare outgoing (or incoming) signaling associated with each ce --------

# Compare outgoing (or incoming) signaling associated with each cell population
# The above analysis summarize the information from the outgoing and incoming signaling together. We can also compare the outgoing (or incoming) signaling pattern between two datasets, allowing to identify signaling pathways/ligand-receptors that exhibit different signaling patterns.

# We can combine all the identified signaling pathways from different datasets and thus compare them side by side, including outgoing signaling, incoming signaling and overall signaling by aggregating outgoing and incoming signaling together. NB: rankNet also shows the comparison of overall signaling, but it does not show the signaling strength in specific cell populations.

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
pathway.union

names(object.list[1])
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
                                         pattern = "outgoing",
                                         color.heatmap = "GnBu",
                                         signaling = pathway.union,
                                         title = names(object.list)[i], width = 6, height = 18)

ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                         pattern = "outgoing",
                                         color.heatmap = "GnBu",
                                         signaling = pathway.union,
                                         title = names(object.list)[i+1], width = 6, height = 18)

pdf("../../out/image/revision/128_comparison_08_signaling_pathway_comparison_outgoing_cellID.pdf",width = 10,height = 10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
dev.off()

#
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
                                         pattern = "incoming",
                                         signaling = pathway.union,
                                         color.heatmap = "GnBu",
                                         title = names(object.list)[i], width = 6, height = 18)

ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                         pattern = "incoming",
                                         color.heatmap = "GnBu",
                                         signaling = pathway.union,
                                         title = names(object.list)[i+1], width = 6, height = 18)

pdf("../../out/image/revision/128_comparison_08_signaling_pathway_comparison_incoming_cellID.pdf",width = 10,height = 10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
dev.off()

# overall
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
                                         pattern = "all",
                                         signaling = pathway.union,
                                         color.heatmap = "GnBu",
                                         title = names(object.list)[i], width = 6, height = 18)

ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                         pattern = "all",
                                         color.heatmap = "GnBu",
                                         signaling = pathway.union,
                                         title = names(object.list)[i+1], width = 6, height = 18)

pdf("../../out/image/revision/128_comparison_08_signaling_pathway_comparison_overall_cellID.pdf",width = 10,height = 10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
dev.off()

# Part III: Identify the upgulated and down-regulated signaling li --------

# Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
# Identify dysfunctional signaling by comparing the communication probabities
# We can compare the communication probabilities mediated by ligand-receptor pairs from some cell groups to other cell groups. This can be done by setting comparison in the function netVisual_bubble.


# -------------------------------------------------------------------------
# # check the position of the levels
# levels(cellchat@idents$joint)
#
# pdf("out/image/09_L_R_pair_IMM_OUT.pdf",width = 8,height = 4)
# # netVisual_bubble(cellchat, sources.use = 1, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
# netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:12), comparison = c(1, 2), angle.x = 45)
# dev.off()
#
# pdf("out/image/09_L_R_pair_ENDO_IN.pdf",width = 8,height = 15)
# # netVisual_bubble(cellchat, sources.use = 1, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
# netVisual_bubble(cellchat, sources.use = c(1:12), targets.use = 6, comparison = c(1, 2), angle.x = 45)
# dev.off()
#
# # why ncam1 ncam2 prob is missing in AC
# subsetCommunication(cellchat, slot.name = "net",
#                     sources.use = c(1:12), targets.use = 6,
#                     signaling = "VEGF")
#
# cellchat.SEN@var.features$features.info %>%
#   dplyr::filter(features %in% c("Vegfa","Kdr"))
#
# cellchat@var.features$CTRL
#
# cellchat@net$CTRL$prob %>%
#   as.data.frame() %>%
#   dplyr::select(contains("Vegfa"))
#
# cellchat@net$CSF$prob %>%
#   as.data.frame() %>%
#   dplyr::select(contains("Vegfa"))
# -------------------------------------------------------------------------

# check the position of the levels
levels(cellchat@idents$joint)

gg1 <- netVisual_bubble(cellchat,
                        sources.use = c("AST","IMM","OLIGO","EXC NEU","INH NEU","VAS","OPC","EPENDYMA","LYM"),
                        targets.use = c("AST","IMM","OLIGO","EXC NEU","INH NEU","VAS","OPC","EPENDYMA","LYM"),
                        comparison = c(1, 2),
                        max.dataset = 2,
                        title.name = "Increased signaling in SEN", angle.x = 45, remove.isolate = T)

gg2 <- netVisual_bubble(cellchat,
                        sources.use = c("AST","IMM","OLIGO","EXC NEU","INH NEU","VAS","OPC","EPENDYMA","LYM"),
                        targets.use = c("AST","IMM","OLIGO","EXC NEU","INH NEU","VAS","OPC","EPENDYMA","LYM"),
                        comparison = c(1, 2),
                        max.dataset = 1,
                        title.name = "Decreased signaling in SEN", angle.x = 45, remove.isolate = T)

pdf("../../out/image/revision/128_comparison_09_L_R_pair_significant_ALL_cellID.pdf",width = 30,height = 15)
gg1 + gg2
dev.off()
#
# NB: The ligand-receptor pairs shown in the bubble plot can be accessed via signaling.LSIncreased = gg1$data.
gg1$data

# Identify dysfunctional signaling by using differential expressio --------

# Identify dysfunctional signaling by using differential expression analysis
# The above method for identifying the upgulated and down-regulated signaling is perfomed by comparing the communication probability between two datasets for each L-R pair and each pair of cell groups. Alternative, we can identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential gene expression analysis. Specifically, we perform differential expression analysis between two biological conditions (i.e., NL and LS) for each cell group, and then obtain the upgulated and down-regulated signaling based on the fCSF change of ligands in the sender cells and receptors in the receiver cells. Such analysis can be done as follows.

# define a positive dataset, i.e., the dataset with positive fCSF change against the other dataset
pos.dataset = "SEN"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat,
                                       group.dataset = "datasets",
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

# net %>% filter(ligand == "Vegfa",receptor == "Kdr")

# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "SEN",ligand.logFC = 0.2, receptor.logFC = NULL)
net.up %>% 
  write_tsv("../../out/table/revision/128_net.up_cellID_SEN.tsv")
# net.up %>%
#   filter(ligand == "Vegfa")
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "CTRL",ligand.logFC = -0.1, receptor.logFC = -0.1)
# net.down2 <- subsetCommunication(cellchat, net = net, datasets = "CTRL",ligand.logFC = -0.1, receptor.logFC = NULL)
# net.down %>%
#   filter(ligand == "Col1a2")

# Since the signaling genes in the net.up and net.down might be complex with multi-subunits, we can do further deconvolution to obtain the individual signaling genes.
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# We then visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.up,
                        sources.use = c("AST","IMM","OLIGO","EXC NEU","INH NEU","VAS","OPC","EPENDYMA","LYM"),
                        targets.use = c("AST","IMM","OLIGO","EXC NEU","INH NEU","VAS","OPC","EPENDYMA","LYM"),
                        comparison = c(1, 2),
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

# gg1 <- netVisual_bubble(cellchat,
#                         pairLR.use = pairLR.use.up,
#                         sources.use = 2,
#                         targets.use = c(1:7),
#                         comparison = c(1, 2),
#                         angle.x = 90,
#                         remove.isolate = T,
#                         title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

# pairLR.use.down = net.down[, "interaction_name", drop = F]
# gg2 <- netVisual_bubble(cellchat,
#                         pairLR.use = pairLR.use.down,
#                         sources.use = c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC"),
#                         targets.use = c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC"),
#                         comparison = c(1, 2),
#                         angle.x = 90,
#                         remove.isolate = T,
#                         title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

pdf("../../out/image/revision/128_comparison_09_L_R_pair_significant_DEG_cellID.pdf",width = 20,height = 8)
# gg1 + gg2
gg1
dev.off()

# Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram

# Chord diagram
levels(cellchat@idents$joint)

# pdf("../../out/image/revision/128_comparison_09_L_R_pair_significant_OLIG_IN_DEG_chord_cellID.pdf",width = 20,height = 20)
# par(mfrow = c(1,1), xpd=TRUE)
# list(netVisual_chord_gene(object.list[[2]],
#                           sources.use = c(1:20),
#                           targets.use = c(11,18),
#                           slot.name = 'net',
#                           net = net.up, lab.cex = 0.8, small.gap = 3.5,
#                           title.name = paste0("Up-regulated signaling in ", names(object.list)[2])),
#      netVisual_chord_gene(object.list[[1]],
#                           sources.use = c(1:20),
#                           targets.use = c(11,18),
#                           slot.name = 'net',
#                           net = net.down, lab.cex = 0.8, small.gap = 3.5,
#                           title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
# )
# dev.off()

# -------------------------------------------------------------------------
# pdf("images/09_L_R_pair_significant_IMM_OUT_DEG_chord.pdf",width = 10,height = 10)
# par(mfrow = c(1,1), xpd=TRUE)
# list(netVisual_chord_gene(object.list[[2]],
#                           pairLR.use = pairLR.use.up,
#                           sources.use = 2,
#                           targets.use = c(1:7),
#                           slot.name = 'net',
#                           # net = net.up,
#                           lab.cex = 0.8, small.gap = 3.5,
#                           title.name = paste0("Up-regulated signaling in ", names(object.list)[2])),
#      netVisual_chord_gene(object.list[[1]],
#                           pairLR.use = pairLR.use.down,
#                           sources.use = 2,
#                           targets.use = c(1:7),
#                           slot.name = 'net',
#                           # net = net.down,
#                           lab.cex = 0.8, small.gap = 3.5,
#                           title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
# )
# dev.off()



# Part IV: Visually compare cell-cell communication using Hierarch --------

# Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
# Similar to the CellChat analysis of individual dataset, we can visualize the cell-cell communication network using Hierarchy plot, Circle plot or Chord diagram.

# Edge color/weight, node color/size/shape: In all visualization plots, edge colors are consistent with the sources as sender, and edge weights are proportional to the interaction strength. Thicker edge line indicates a stronger signal. In the Hierarchy plot and Circle plot, circle sizes are proportional to the number of cells in each cell group. In the hierarchy plot, solid and open circles represent source and target, respectively. In the Chord diagram, the inner thinner bar colors represent the targets that receive signal from the corresponding outer bar. The inner bar size is proportional to the signal strength received by the targets. Such inner bar is helpful for interpreting the complex chord diagram. Note that there exist some inner bars without any chord for some cell groups, please just igore it because this is an issue that has not been addressed by circlize package.

# # choose a pathways
# sort(table(net.up %>%
#              filter(source =="IMM"|target =="AST") %>%
#              pull(pathway_name)),decreasing = T)
# sort(table(net.down %>%
#              filter(source =="IMM"|target =="AST") %>%
#              pull(pathway_name)),decreasing = T)
#
# pathways.show <- c("LAMININ")
#
# # control the edge weights across different datasets
# weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
#
# pdf("images/09_circle_LAMININ.pdf",width = 10,height = 5)
# par(mfrow = c(1,2), xpd=TRUE)
# lapply(1:length(object.list),function(i){
#   netVisual_aggregate(object.list[[i]],
#                       signaling = pathways.show,
#                       layout = "circle",
#                       edge.weight.max = weight.max[1],
#                       edge.width.max = 10,
#                       signaling.name = paste(pathways.show, names(object.list)[i]))
# })
# dev.off()
#
# #
# pathways.show <- c("PTN")
#
# # control the edge weights across different datasets
# weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
#
# pdf("images/09_circle_PTN.pdf",width = 10,height = 5)
# par(mfrow = c(1,2), xpd=TRUE)
# lapply(1:length(object.list),function(i){
#   netVisual_aggregate(object.list[[i]],
#                       signaling = pathways.show,
#                       layout = "circle",
#                       edge.weight.max = weight.max[1],
#                       edge.width.max = 10,
#                       signaling.name = paste(pathways.show, names(object.list)[i]))
# })
# dev.off()
#
# # -------------------------------------------------------------------------
# pathways.show <- c("COLLAGEN")
# par(mfrow = c(1,2), xpd=TRUE)
#
# ht <- lapply(1:length(object.list),function(i){
#   netVisual_heatmap(object.list[[i]],
#                     signaling = pathways.show,
#                     color.heatmap = "Reds",
#                     title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
#
# })
#
# pdf("images/09_heatmap_COLLAGEN.pdf",width = 10,height = 5)
# ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# dev.off()
#
# pathways.show <- c("PTN")
# par(mfrow = c(1,2), xpd=TRUE)
#
# ht <- lapply(1:length(object.list),function(i){
#   netVisual_heatmap(object.list[[i]],
#                     signaling = pathways.show,
#                     color.heatmap = "Reds",
#                     title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
#
# })
#
# pdf("images/09_heatmap_PTN.pdf",width = 10,height = 5)
# ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# dev.off()
#
# # -------------------------------------------------------------------------
# # Chord diagram
# pathways.show <- c("COLLAGEN")
#
# pdf("images/09_chord_COLLAGEN.pdf",width = 10,height = 5)
# par(mfrow = c(1,2), xpd=TRUE)
# lapply(1:length(object.list),function(i){
#   netVisual_aggregate(object.list[[i]],
#                       signaling = pathways.show,
#                       layout = "chord",
#                       signaling.name = paste(pathways.show, names(object.list)[i]))
# })
# dev.off()
#
# pathways.show <- c("PTN")
#
# pdf("images/09_chord_PTN.pdf",width = 10,height = 5)
# par(mfrow = c(1,2), xpd=TRUE)
# lapply(1:length(object.list),function(i){
#   netVisual_aggregate(object.list[[i]],
#                       signaling = pathways.show,
#                       layout = "chord",
#                       signaling.name = paste(pathways.show, names(object.list)[i]))
# })
# dev.off()

# For the chord diagram, CellChat has an independent function netVisual_chord_cell to flexibly visualize the signaling network by adjusting different parameters in the circlize package. For example, we can define a named char vector group to create multiple-group chord diagram, e.g., grouping cell clusters into different cell types.

# Chord diagram
# levels(cellchat@idents$joint)
# grouping cell clusters into fibroblast, DC and TC cells
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
# names(group.cellType) <- levels(object.list[[1]]@idents)
# pathways.show <- c("CXCL")
#
# pdf("09_chord_grouped.pdf",width = 10,height = 5)
# par(mfrow = c(1,2), xpd=TRUE)
# lapply(1:length(object.list),function(i){
#   netVisual_chord_cell(object.list[[i]],
#                        signaling = pathways.show,
#                        group = group.cellType,
#                        title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
# })
# dev.off()

# Using chord diagram, CellChat provides two functions netVisual_chord_cell and netVisual_chord_gene for visualizing cell-cell communication with different purposes and different levels. netVisual_chord_cell is used for visualizing the cell-cell communication between different cell groups (where each sector in the chord diagram is a cell group), and netVisual_chord_gene is used for visualizing the cell-cell communication mediated by mutiple ligand-receptors or signaling pathways (where each sector in the chord diagram is a ligand, receptor or signaling pathway.)

# # compare all the interactions sending from Inflam.FIB to DC cells
# pdf("09_L_R_pair_significant_DEG_chord_subset.pdf",width = 10,height = 10)
# par(mfrow = c(1,1), xpd=TRUE)
# lapply(1:length(object.list),function(i){
#   netVisual_chord_gene(object.list[[i]],
#                        sources.use = 4,
#                        targets.use = c(5:8),
#                        lab.cex = 0.5,
#                        title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
# })
# dev.off()

# compare all the interactions sending from fibroblast to inflamatory immune cells
# pdf("09_L_R_pair_significant_DEG_chord_subset2.pdf",width = 10,height = 10)
# par(mfrow = c(1,1), xpd=TRUE)
# lapply(1:length(object.list),function(i){
#   netVisual_chord_gene(object.list[[i]],
#                        sources.use = c(1,2, 3, 4),
#                        targets.use = c(8,10),
#                        title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
# })
# dev.off()

# show all the significant signaling pathways from fibroblast to immune cells
# pdf("09_L_R_pair_significant_DEG_chord_subset3.pdf",width = 10,height = 10)
# par(mfrow = c(1,1), xpd=TRUE)
# lapply(1:length(object.list),function(i){
#   netVisual_chord_gene(object.list[[i]],
#                        sources.use = c(1,2,3,4),
#                        targets.use = c(5:11),
#                        slot.name = "netP",
#                        title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
# })
# dev.off()

# NB: Please ignore the note when generating the plot such as “Note: The first link end is drawn out of sector ‘MIF’.”. If the gene names are  overlapped, you can adjust the argument small.gap by decreasing the value.


#  Part V: Compare the signaling gene expression distribution betw --------
#  Part V: Compare the signaling gene expression distribution between different datasets
# We can plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.

# set factor level
# net.up %>%
# net.down %>%

# -------------------------------------------------------------------------
# net %>%
#   # filter(str_detect(pathway_name,pattern = "FN"))
#   # filter(str_detect(interaction_name,pattern = "FN"))
#   filter(ligand == "Tgfa")
#
# cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = c("CTRL", "CSF"))
#
# pdf("out/image/10_violin_Tgfa.pdf",width = 10,height = 10)
# plotGeneExpression(cellchat, signaling = "EGF", split.by = "datasets", colors.ggplot = T)
# dev.off()

# -------------------------------------------------------------------------
# focus on the APP signalling

# list(netVisual_chord_gene(object.list[[2]],
#                           sources.use = 5,
#                           targets.use = 1:7,
#                           slot.name = 'netP',
#                           net = net.up, lab.cex = 0.8, small.gap = 3.5,
#                           title.name = paste0("Up-regulated signaling in ", names(object.list)[2])),
#      netVisual_chord_gene(object.list[[1]],
#                           sources.use = 5,
#                           targets.use = 1:7,
#                           slot.name = 'netP',
#                           net = net.down, lab.cex = 0.8, small.gap = 3.5,
#                           title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
# )

# # -------------------------------------------------------------------------
# subsetCommunication(cellchat,signaling = "EGF")
# # return the dataset of DE for the communication
# # test <- identifyOverExpressedGenes(cellchat,
# #                                    group.dataset = "datasets",
# #                                    pos.dataset = pos.dataset,
# #                                    features.name = features.name,
# #                                    only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1,return.object = F)
#
# net <- netMappingDEG(cellchat, features.name = features.name)
# net.up <- subsetCommunication(cellchat, net = net, datasets = "CSF",ligand.logFC = 0.2, receptor.logFC = NULL)
#
# net %>%
#   filter(ligand == "Tgfa")
#
# # net.up2 %>%
# #   filter(ligand == "APP")
#
# list(
#   netVisual_aggregate(object.list[[2]], signaling = "EGF", layout = "circle", signaling.name = paste("EGF", names(object.list)[2])),
#   netVisual_aggregate(object.list[[1]], signaling = "EGF", layout = "circle", signaling.name = paste("EGF", names(object.list)[1]))
# )
#
# netVisual_chord_gene(object.list[[2]],
#                      slot.name = 'net',
#                      net = net.up, lab.cex = 0.8, small.gap = 3.5,
#                      title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
# # -------------------------------------------------------------------------

# # here is  a way to summarise the numkber of LR pairs significantly disregulated up and down
# # netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#
# groupSize <- as.numeric(table(cellchat@idents$joint))
#
# pdf("out/image/00_test.pdf",width = 7,height = 7)
# netVisual_circle(net = table(net.up$source,net.up$target),vertex.weight = groupSize,weight.scale = T,
#                  label.edge= F, title.name = "DE UP Number of interactions")
# dev.off()
#
# pdf("out/image/00_test_red.pdf",width = 7,height = 7)
# netVisual_circle2(net = table(net.up$source,net.up$target),vertex.weight = groupSize,weight.scale = T,
#                   label.edge= F, title.name = "DE UP Number of interactions")
# dev.off()

# -------------------------------------------------------------------------
# APP CD74
# plot the network for the specific LR pair in eachc individula dataset
# confirm it by checking alos the expression table
# check the level of expression per cell type and dataset
test <- paste0(cellchat@meta$expertAnno.l1,"_",cellchat@meta$datasets)
# add it to the metadata
cellchat@meta$test <- test

unique(c(unique(net.up$ligand),unique(net.up$receptor)))

computeAveExpr(cellchat, features = c("CD44","SPP1"),group.by = c("test"))
mat_test <- computeAveExpr(cellchat, features = unique(c(unique(net.up$ligand),unique(net.up$receptor))),group.by = c("test"))

net.up

mat_test_general <- computeAveExpr(cellchat, features = unique(c(unique(net.up$ligand),unique(net.up$receptor))),group.by = c("expertAnno.l1")) %>%
  .[,c("AST","EPENDYMA","EXC NEU","INH NEU","LYM","OLIGO","OPC","VAS")]

# add the avg expression for imm senescent and control
mat_test_update <- left_join(mat_test %>% data.frame() %>% select("IMM_CTRL","IMM_SEN") %>% rownames_to_column("gene"),
                             mat_test_general %>% data.frame() %>% rownames_to_column("gene"),by = "gene") %>%
  column_to_rownames("gene")


pdf("../../out/image/revision/128_comparison_09_L_R_pair_significant_DEG_cellID_heatmap_avgExp.pdf",width = 6,height = 7)
Heatmap(mat_test_update,viridis::viridis(option = "turbo",n = 10),cluster_columns = F,name = "avg_exp")
dev.off()

# try to scale by gene
t(scale(t(mat_test)))
rowSums(t(scale(t(mat_test))))

Heatmap(t(scale(t(mat_test))),
        # viridis::viridis(option = "turbo",n = 10),
        col = colorRamp2(c(-2, 0, 3), c("blue", "white", "red")),
        cluster_columns = F,name = "scaled_exp")

# separate ligand and receptors
mat_test_ligand <- computeAveExpr(cellchat, features = unique(net.up$ligand),group.by = c("test"))
mat_test_receptor <- computeAveExpr(cellchat, features = unique(net.up$receptor),group.by = c("test"))


hm_ligand <- Heatmap(t(scale(t(mat_test_ligand))),
                     # viridis::viridis(option = "turbo",n = 10),
                     col = colorRamp2(c(-2, 0, 3), c("blue", "white", "red")),
                     cluster_columns = F,name = "scaled_exp",row_title = "ligand")

hm_receptor <- Heatmap(t(scale(t(mat_test_receptor))),
                       # viridis::viridis(option = "turbo",n = 10),
                       col = colorRamp2(c(-2, 0, 3), c("blue", "white", "red")),
                       cluster_columns = F,name = "scaled_exp",row_title = "receptor")

hm_ligand %v% hm_receptor

# Martina suggested to generate the matrices for the other cell types not split by senescence or control. remove the IMM.
mat_test_ligand_general <- computeAveExpr(cellchat, features = unique(net.up$ligand),group.by = c("expertAnno.l1")) %>%
  .[,c("AST","EPENDYMA","EXC NEU","INH NEU","LYM","OLIGO","OPC","VAS")]
mat_test_receptor_general <- computeAveExpr(cellchat, features = unique(net.up$receptor),group.by = c("expertAnno.l1")) %>%
  .[,c("AST","EPENDYMA","EXC NEU","INH NEU","LYM","OLIGO","OPC","VAS")]

# add the avg expression for imm senescent and control
mat_test_ligand_update <- left_join(mat_test_ligand %>% data.frame() %>% select("IMM_CTRL","IMM_SEN") %>% rownames_to_column("gene"),
          mat_test_ligand_general %>% data.frame() %>% rownames_to_column("gene"),by = "gene") %>%
  column_to_rownames("gene")

mat_test_receptor_update <- left_join(mat_test_receptor %>% data.frame() %>% select("IMM_CTRL","IMM_SEN") %>% rownames_to_column("gene"),
                                              mat_test_receptor_general %>% data.frame() %>% rownames_to_column("gene"),by = "gene") %>%
  column_to_rownames("gene")

hm_ligand_update <- Heatmap(t(scale(t(mat_test_ligand_update))),
                            # viridis::viridis(option = "turbo",n = 10),
                            col = colorRamp2(c(-2, 0, 3), c("blue", "white", "red")),
                            cluster_columns = F,name = "scaled_exp",row_title = "ligand")

hm_receptor_update <- Heatmap(t(scale(t(mat_test_receptor_update))),
                              # viridis::viridis(option = "turbo",n = 10),
                              col = colorRamp2(c(-2, 0, 3), c("blue", "white", "red")),
                              cluster_columns = F,name = "scaled_exp",row_title = "receptor")

pdf("../../out/image/revision/128_comparison_09_L_R_pair_significant_DEG_cellID_heatmap_scaledExp.pdf",width = 6,height = 7)
hm_ligand_update %v% hm_receptor_update
dev.off()

# pdf("out/image/10_LR_pairs_APP.pdf",width = 7,height = 7)
# list(
#   netVisual_chord_gene2(object.list[[2]],pairLR.use = data.frame('interaction_name' = "BMP5_ACVR1_BMPR2"),
#                        slot.name = 'net',
#                        net = net %>%
#                          filter(datasets=="CSF"), lab.cex = 0.8, small.gap = 3.5,
#                        title.name = paste0("significant signaling in ", names(object.list)[2])),
#
#   netVisual_chord_gene2(object.list[[1]],pairLR.use = data.frame('interaction_name' = "BMP5_ACVR1_BMPR2"),
#                        slot.name = 'net',
#                        net = net %>%
#                          filter(datasets=="CTRL"), lab.cex = 0.8, small.gap = 3.5,
#                        title.name = paste0("significant signaling in ", names(object.list)[1]))
# )
# dev.off()
# plot only the connection that are upregulated in CCA
# pdf("images/10_LR_pairs_APP_UP_CCA.pdf",width = 7,height = 7)
# list(
#   netVisual_chord_gene2(object.list[[2]],pairLR.use = data.frame('interaction_name' = "APP_CD74"),
#                        slot.name = 'net',
#                        net = net.up, lab.cex = 0.8, small.gap = 3.5,
#                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
# )
# dev.off()
# -------------------------------------------------------------------------
# whole dataset
# plot all the pathway that are goind up in CCA
pdf("../../out/image/revision/128_comparison_09_pathway_significant_CSF_chord_cellID.pdf",width = 15,height = 15)
list(
  netVisual_chord_gene(object.list[[2]],
                        slot.name = 'netP',
                        net = net.up, lab.cex = 0.8, small.gap = 3.5,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
  # netVisual_chord_gene(object.list[[1]],
  #                       slot.name = 'netP',
  #                       net = net.down, lab.cex = 0.8, small.gap = 3.5,
  #                       title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
)
dev.off()

# check the position of the levels
levels(cellchat@idents$joint)

# netVisual_chord_gene(object.list[[2]],
#                      sources.use = c(1:20),
#                      targets.use = c(11,18),
#                      slot.name = 'net',
#                      net = net.up, lab.cex = 0.8, small.gap = 3.5,
#                      title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

# pdf("../../out/image/revision/128_comparison_09_pathway_significant_CSF_chord_OLIG_IN_cellID.pdf",width = 15,height = 15)
# list(
#   netVisual_chord_gene2(object.list[[2]],targets.use = c(17),
#                         slot.name = 'netP',
#                         net = net.up, lab.cex = 0.8, small.gap = 3.5,
#                         title.name = paste0("Up-regulated signaling in ", names(object.list)[2])),
#   netVisual_chord_gene2(object.list[[1]],targets.use = c(17),
#                         slot.name = 'netP',
#                         net = net.down, lab.cex = 0.8, small.gap = 3.5,
#                         title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
# )
# dev.off()

# pdf("out/image/09_pathway_significant_CSF_chord_OLIG_IN2.pdf",width = 15,height = 15)
# list(
#   netVisual_chord_gene2(object.list[[2]],targets.use = c(11,18),
#                         slot.name = 'net',
#                         net = net.up, lab.cex = 0.8, small.gap = 3.5,
#                         title.name = paste0("Up-regulated signaling in ", names(object.list)[2])),
#   netVisual_chord_gene2(object.list[[1]],targets.use = c(11,18),
#                         slot.name = 'net',
#                         net = net.down, lab.cex = 0.8, small.gap = 3.5,
#                         title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
# )
# dev.off()


# plot all the LR
pdf("../../out/image/revision/128_comparison_09_LR_significant_CSF_chord_cellID.pdf",width = 15,height = 15)
list(
  netVisual_chord_gene(object.list[[2]],
                        slot.name = 'net',
                        net = net.up, lab.cex = 0.8, small.gap = 3.5,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
  # netVisual_chord_gene2(object.list[[1]],
  #                       slot.name = 'net',
  #                       net = net.down, lab.cex = 0.8, small.gap = 3.5,
  #                       title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
)
dev.off()

netVisual_chord_gene(object.list[[2]],
                     slot.name = 'net',
                     net = net.up, lab.cex = 1, small.gap = 10,big.gap = 100, scale = T,
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

pdf("../../out/image/revision/128_comparison_09_LR_significant_CSF_chord_cellID_high.pdf",width = 10,height = 10)
list(
  netVisual_chord_gene(object.list[[2]],
                       slot.name = 'net',
                       net = net.up, lab.cex = 0.7, small.gap = 0.01,
                       title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
  # netVisual_chord_gene2(object.list[[1]],
  #                       slot.name = 'net',
  #                       net = net.down, lab.cex = 0.8, small.gap = 3.5,
  #                       title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
)
dev.off()

pdf("../../out/image/revision/128_comparison_09_LR_significant_CSF_chord_cellID_high2.pdf",width = 7,height = 7)
list(
  netVisual_chord_gene(object.list[[2]],
                       slot.name = 'net',
                       net = net.up, lab.cex = 0.8, small.gap = 1,big.gap = 1,scale = T,
                       title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
  # netVisual_chord_gene2(object.list[[1]],
  #                       slot.name = 'net',
  #                       net = net.down, lab.cex = 0.8, small.gap = 3.5,
  #                       title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
)
dev.off()

# Save the merged CellChat object -----------------------------------------
saveRDS(cellchat, file = "../../out/object/cellchat_comparison_CTRL_vs_CSF_cellID.rds")
