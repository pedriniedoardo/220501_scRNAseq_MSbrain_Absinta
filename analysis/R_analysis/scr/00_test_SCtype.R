# # es.max <- sctype_score(scRNAseqData = scobj[["RNA"]]@scale.data, scaled = TRUE, 
# #                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# 
# sctype_score
# 
# scRNAseqData = scobj[["RNA"]]@scale.data
# scRNAseqData[1:10,1:10]
# scRNAseqData1 = scobj[["RNA"]]@data
# scRNAseqData1[1:10,1:10]
# 
# scaled = F, 
# gs = gs_list$gs_positive
# gs2 = gs_list$gs_negative
# 
# # -------------------------------------------------------------------------
# # marker sensitivity
# marker_stat = sort(table(unlist(gs)), decreasing = T); 
# marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
#                                 gene_ = names(marker_stat), stringsAsFactors = !1)
# 
# # # convert gene names to Uppercase
# # if(gene_names_to_uppercase){
# #   rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
# # }
# 
# # subselect genes only found in data
# names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
# 
# rownames(scRNAseqData)
# 
# 
# gs = lapply(1:length(gs), function(d_){ 
#   GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
# 
# # lapply(names_gs_cp, function(d_){
# #   GeneIndToKeep = rownames(scRNAseqData) %in% unlist(gs[[d_]])
# #   rownames(scRNAseqData)[GeneIndToKeep]
# #   }) %>% 
# #   setNames(names_gs_cp)
# 
# gs2 = lapply(1:length(gs2), function(d_){ 
#   GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
# 
# names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
# 
# cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
# 
# # z-scale if not
# if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
# 
# # multiple by marker sensitivity
# for(jj in 1:nrow(cell_markers_genes_score)){
#   Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
# }
# 
# # subselect only with marker genes
# Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
# 
# # combine scores
# es = do.call("rbind", lapply(names(gs), function(gss_){ 
#   sapply(1:ncol(Z), function(j) {
#     gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
#     sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
#     if(is.na(sum_t2)){
#       sum_t2 = 0;
#     }
#     sum_t1 + sum_t2
#   })
# })) 
# 
# dimnames(es) = list(names(gs), colnames(Z))
# es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
# 
# es.max