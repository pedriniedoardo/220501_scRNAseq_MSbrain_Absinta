# AIM ---------------------------------------------------------------------
# follow up on scr/40_plotGene_integrationSkip_manualClean_harmony_update.R
# TSPO: test the pos vs neg TSPO classification of MG (IMM) and VAS cells danny's concern: the pos/neg TSPO split could simply be tracking cell quality/sequencing depth rather than real biology
# 1) test whether nCount_RNA/nFeature_RNA distributions are homogeneous between TSPO pos and TSPO neg cells, for MG (IMM) and VAS
# 2) correlate continuous TSPO expression with sequencing depth, at the cell level and at the sample (pseudobulk) level
# 3) quantify with a logistic regression how well nCount_RNA/nFeature_RNA alone predict the pos/neg classification
# 4) build a depth-matched pos/neg comparison (nearest-neighbour matching on nCount_RNA/nFeature_RNA) and re-run the DGE on the matched cells, which directly tests whether the classification itself (not just the downstream fold-change estimate) is confounded by depth
# 5) repeat the TSPO pos vs neg DGE splitting the cells by disease status (CTRL vs MS, metadata column "disease")

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(RANN)
library(broom)
library(pROC)
library(ggridges)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")

GOI <- "TSPO"

table(data.combined$disease)

# helper: subset a cell type, tag TSPO pos/neg -----------------------------
add_TSPO_cat <- function(sobj){
  exp_TSPO <- FetchData(sobj,GOI,slot = "data") %>%
    rownames_to_column("barcodes")
  
  sobj$TSPO <- sobj@meta.data %>%
    rownames_to_column("barcodes") %>%
    left_join(exp_TSPO,by = "barcodes") %>%
    pull(GOI)
  
  sobj$TSPO_cat <- case_when(sobj$TSPO == 0~"neg",
                             T~"pos")
  sobj
}

DefaultAssay(data.combined) <- "RNA"

sobj_MG <- subset(data.combined,subset = expertAnno.l1 == 'IMM') %>% add_TSPO_cat()
sobj_VAS <- subset(data.combined,subset = expertAnno.l1 == 'VAS') %>% add_TSPO_cat()

table(sobj_MG$TSPO_cat)
table(sobj_VAS$TSPO_cat)

list_sobj <- list(MG = sobj_MG, VAS = sobj_VAS)

# TSPO QC homogeneity test: nCount_RNA / nFeature_RNA pos vs neg -----------
df_QC_TSPO <- imap(list_sobj,function(sobj,cell_id){
  sobj@meta.data %>%
    dplyr::select(nCount_RNA,nFeature_RNA,TSPO_cat) %>%
    mutate(cell_id = cell_id)
}) %>%
  bind_rows()

# save the QC table
df_QC_TSPO %>%
  write_tsv("../../out/table/custom/07_TSPO_df_QC_posVSneg.tsv")

# visualize the distributions
df_QC_TSPO %>%
  pivot_longer(names_to = "metric",values_to = "value",cols = c(nCount_RNA,nFeature_RNA)) %>%
  ggplot(aes(y = TSPO_cat,x = value)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha=0.5, vline_linetype="dashed") +
  # geom_violin(scale = "width") +
  # geom_boxplot(width = 0.1,outlier.shape = NA) +
  facet_wrap(cell_id~metric,scales = "free_x") +
  scale_x_log10() +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("../../out/image/custom/07_TSPO_violin_QC_posVSneg.pdf",width = 7,height = 9)

# test homogeneity of the two distributions per metric and cell type wilcox.test: tests whether the two groups have the same location (median)
# ks.test: tests whether the two groups come from the same distribution (shape+location)
df_QC_TSPO_test <- df_QC_TSPO %>%
  group_by(cell_id) %>%
  summarise(
    nCount_wilcox_p = wilcox.test(nCount_RNA~TSPO_cat)$p.value,
    nCount_ks_p = ks.test(nCount_RNA[TSPO_cat=="pos"],nCount_RNA[TSPO_cat=="neg"])$p.value,
    nFeature_wilcox_p = wilcox.test(nFeature_RNA~TSPO_cat)$p.value,
    nFeature_ks_p = ks.test(nFeature_RNA[TSPO_cat=="pos"],nFeature_RNA[TSPO_cat=="neg"])$p.value,
    .groups = "drop")

df_QC_TSPO_test %>%
  write_tsv("../../out/table/custom/07_TSPO_df_QC_test_posVSneg.tsv")

df_QC_TSPO_test

# TSPO expression vs sequencing depth ---------------------------------------
# point 1: does TSPO detection just track library complexity?
# cell-level: correlate the normalized TSPO expression with nCount/nFeature
df_TSPO_depth <- imap(list_sobj,function(sobj,cell_id){
  sobj@meta.data %>%
    rownames_to_column("barcodes") %>%
    dplyr::select(barcodes,orig.ident,nCount_RNA,nFeature_RNA,TSPO,TSPO_cat) %>%
    mutate(cell_id = cell_id)
}) %>%
  bind_rows()

df_TSPO_depth_cor <- df_TSPO_depth %>%
  group_by(cell_id) %>%
  summarise(
    cor_TSPO_nCount = cor(TSPO,nCount_RNA,method = "spearman"),
    cor_TSPO_nFeature = cor(TSPO,nFeature_RNA,method = "spearman"),
    .groups = "drop")

df_TSPO_depth_cor %>%
  write_tsv("../../out/table/custom/07_TSPO_correlation_expression_depth_cellLevel.tsv")

df_TSPO_depth_cor

df_TSPO_depth %>%
  ggplot(aes(x = nCount_RNA,y = TSPO)) +
  geom_point(alpha = 0.1,size = 0.3) +
  geom_smooth(method = "lm",col = "red") +
  facet_wrap(~cell_id,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/custom/07_TSPO_scatter_expression_vs_nCount_cellLevel.pdf",width = 8,height = 4)

# sample level (pseudobulk): does the proportion of TSPO+ cells per sample
# simply track that sample's sequencing depth?
df_TSPO_sample <- df_TSPO_depth %>%
  group_by(cell_id,orig.ident) %>%
  summarise(n_cells = n(),
            prop_pos = mean(TSPO_cat == "pos"),
            median_nCount = median(nCount_RNA),
            median_nFeature = median(nFeature_RNA),
            .groups = "drop")

df_TSPO_sample %>%
  write_tsv("../../out/table/custom/07_TSPO_df_sample_propPos_vs_depth.tsv")

df_TSPO_sample %>%
  group_by(cell_id) %>%
  summarise(cor_propPos_nCount = cor(prop_pos,median_nCount,method = "spearman"),
            cor_propPos_nFeature = cor(prop_pos,median_nFeature,method = "spearman"),
            .groups = "drop") %>%
  write_tsv("../../out/table/custom/07_TSPO_correlation_propPos_depth_sampleLevel.tsv")

df_TSPO_sample %>%
  ggplot(aes(x = median_nCount,y = prop_pos)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = orig.ident),size = 2) +
  geom_smooth(method = "lm",col = "red") +
  facet_wrap(~cell_id,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/custom/07_TSPO_scatter_sample_propPos_vs_depth.pdf",width = 8,height = 4)

# logistic regression: how well does depth alone predict TSPO_cat -----------
# point 3: quantify how much of the pos/neg split is explained by depth alone
df_TSPO_logit <- imap(list_sobj,function(sobj,cell_id){
  glm_dat <- sobj@meta.data %>%
    mutate(TSPO_pos = as.integer(TSPO_cat == "pos"))
  
  fit <- glm(TSPO_pos~nCount_RNA+nFeature_RNA,data = glm_dat,family = binomial)
  fit_null <- glm(TSPO_pos~1,data = glm_dat,family = binomial)
  
  pseudo_R2 <- as.numeric(1 - (logLik(fit)/logLik(fit_null)))
  model_auc <- as.numeric(auc(roc(glm_dat$TSPO_pos,fitted(fit),quiet = TRUE)))
  
  broom::tidy(fit) %>%
    mutate(cell_id = cell_id,
           pseudo_R2 = pseudo_R2,
           AUC = model_auc)
}) %>%
  bind_rows()

df_TSPO_logit %>%
  write_tsv("../../out/table/custom/07_TSPO_logit_classification_by_depth.tsv")

df_TSPO_logit

# depth-matched pos vs neg comparison ----------------------------------------
# point 5: nearest-neighbour match TSPO+ and TSPO- cells on nCount_RNA/nFeature_RNA so the two groups are depth-balanced by construction, then re-check homogeneity and re-run the DGE on the matched subset. this tests whether the classification itself (not just the fold-change estimate) is confounded by depth, which adjusting for nCount_RNA as a DE latent variable would not address.
# sobj <- list_sobj$VAS
# k_search <- 50
match_depth <- function(sobj,k_search = 50){
  meta <- sobj@meta.data %>%
    rownames_to_column("barcodes")
  
  # scale on log1p to stabilize the variance of the count-based covariates
  X <- meta %>%
    transmute(log_nCount = log1p(nCount_RNA),
              log_nFeature = log1p(nFeature_RNA)) %>%
    scale()
  
  grp <- meta$TSPO_cat
  tab <- table(grp)
  # match to the size of the smaller group
  query_lab <- names(tab)[which.min(tab)]
  pool_lab <- names(tab)[which.max(tab)]
  
  query_idx <- which(grp == query_lab)
  pool_idx <- which(grp == pool_lab)
  
  k_use <- min(k_search,length(pool_idx))
  nn_res <- nn2(data = X[pool_idx,,drop = FALSE],query = X[query_idx,,drop = FALSE],k = k_use)
  
  # greedy assignment: resolve the best-matched query cells first, without replacement on the pool side
  used <- rep(FALSE,length(pool_idx))
  matched_pool_local <- rep(NA_integer_,length(query_idx))
  order_q <- order(nn_res$nn.dists[,1])
  
  for(i in order_q){
    for(j in seq_len(k_use)){
      cand <- nn_res$nn.idx[i,j]
      if(!used[cand]){
        used[cand] <- TRUE
        matched_pool_local[i] <- cand
        break
      }
    }
  }
  
  keep <- !is.na(matched_pool_local)
  matched_barcodes <- c(meta$barcodes[query_idx[keep]],
                        meta$barcodes[pool_idx[matched_pool_local[keep]]])
  
  subset(sobj,cells = matched_barcodes)
}

list_sobj_matched <- map(list_sobj,match_depth)

table(list_sobj_matched$MG$TSPO_cat)
table(list_sobj_matched$VAS$TSPO_cat)

# confirm the matched groups are now depth-balanced
df_QC_TSPO_matched <- imap(list_sobj_matched,function(sobj,cell_id){
  sobj@meta.data %>%
    dplyr::select(nCount_RNA,nFeature_RNA,TSPO_cat) %>%
    mutate(cell_id = cell_id)
}) %>%
  bind_rows()

df_QC_TSPO_matched_test <- df_QC_TSPO_matched %>%
  group_by(cell_id) %>%
  summarise(
    nCount_wilcox_p = wilcox.test(nCount_RNA~TSPO_cat)$p.value,
    nCount_ks_p = ks.test(nCount_RNA[TSPO_cat=="pos"],nCount_RNA[TSPO_cat=="neg"])$p.value,
    nFeature_wilcox_p = wilcox.test(nFeature_RNA~TSPO_cat)$p.value,
    nFeature_ks_p = ks.test(nFeature_RNA[TSPO_cat=="pos"],nFeature_RNA[TSPO_cat=="neg"])$p.value,
    .groups = "drop"
  )

df_QC_TSPO_matched_test %>%
  write_tsv("../../out/table/custom/07_TSPO_df_QC_test_posVSneg_depthMatched.tsv")

df_QC_TSPO_matched_test

df_QC_TSPO_matched %>%
  pivot_longer(names_to = "metric",values_to = "value",cols = c(nCount_RNA,nFeature_RNA)) %>%
  ggplot(aes(y = TSPO_cat,x = value)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha=0.5, vline_linetype="dashed") +
  # geom_violin(scale = "width") +
  # geom_boxplot(width = 0.1,outlier.shape = NA) +
  facet_wrap(cell_id~metric,scales = "free_x") +
  scale_x_log10() +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("../../out/image/custom/07_TSPO_violin_QC_posVSneg_depthMatched.pdf",width = 9,height = 7)

# re-run the pos vs neg DGE on the depth-matched cells
run_DGE_TSPO <- function(sobj,cell_id){
  Idents(sobj) <- "TSPO_cat"
  RunPresto(object = sobj,ident.1 = "pos",ident.2 = "neg",min.pct = 0.01) %>%
    rownames_to_column("gene") %>%
    mutate(cell_id = cell_id)
}

df_DGE_TSPO_matched <- imap(list_sobj_matched,run_DGE_TSPO) %>%
  bind_rows()

df_DGE_TSPO_matched %>%
  write_tsv("../../out/table/custom/07_TSPO_res_test_posVSneg_depthMatched.tsv")

# compare the matched vs the original (unmatched) full DGE, per cell type, to see whether the DE signal survives depth matching
df_DGE_TSPO_full <- imap(list_sobj,run_DGE_TSPO) %>%
  bind_rows()

df_DGE_compare <- df_DGE_TSPO_full %>%
  dplyr::select(cell_id,gene,avg_log2FC_full = avg_log2FC,p_val_adj_full = p_val_adj) %>%
  inner_join(df_DGE_TSPO_matched %>%
               dplyr::select(cell_id,gene,avg_log2FC_matched = avg_log2FC,p_val_adj_matched = p_val_adj),
             by = c("cell_id","gene")) %>%
  mutate(DE_cat_full = case_when(avg_log2FC_full > 0.5 & p_val_adj_full < 0.01~"up",
                                 avg_log2FC_full < (-0.5) & p_val_adj_full < 0.01~"down",
                                 T~"no")) %>%
  mutate(DE_cat_matched = case_when(avg_log2FC_matched > 0.5 & p_val_adj_matched < 0.01~"up",
                                    avg_log2FC_matched < (-0.5) & p_val_adj_matched < 0.01~"down",
                                    T~"no"))

df_DGE_compare %>%
  write_tsv("../../out/table/custom/07_TSPO_DGE_compare_full_vs_depthMatched.tsv")

df_DGE_compare %>%
  group_by(cell_id) %>%
  summarise(cor_log2FC = cor(avg_log2FC_full,avg_log2FC_matched,method = "spearman"),
            .groups = "drop")

# add confusion matrix for significance assignament
df_DGE_compare %>%
  group_by(cell_id,DE_cat_full,DE_cat_matched) %>%
  summarise(n = n(),.groups = "drop") %>%
  group_by(cell_id) %>%
  mutate(tot = sum(n),
         prop = n/tot) %>%
  arrange(cell_id,DE_cat_full,DE_cat_matched)

# shwo the comparisoin anly fro the one significant in either of the two
df_DGE_compare %>%
  filter(!(DE_cat_full == "no" & DE_cat_matched == "no")) %>%
  ggplot(aes(x = avg_log2FC_full,y = avg_log2FC_matched)) +
  geom_point(alpha = 0.3,size = 0.5) +
  geom_abline(slope = 1,intercept = 0,linetype = "dashed",col = "red") +
  facet_wrap(~cell_id,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/custom/07_TSPO_scatter_log2FC_full_vs_depthMatched.pdf",width = 8,height = 4)

df_DGE_compare %>%
  filter(!(DE_cat_full == "no" & DE_cat_matched == "no")) %>%
  ggplot(aes(x = p_val_adj_full,y = p_val_adj_matched)) +
  geom_point(alpha = 0.3,size = 0.5) +
  geom_hline(yintercept = 0.05,linetype = "dashed",col = "red") +
  geom_vline(xintercept = 0.05,linetype = "dashed",col = "red") +
  facet_wrap(~cell_id,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/custom/07_TSPO_scatter_padj_full_vs_depthMatched.pdf",width = 8,height = 4)

# random (depth-agnostic) same-size null for the matched comparison --------
# the depth-matched groups are smaller than the full dataset, so comparing significance calls directly confounds "lost power from smaller n" with "lost signal from removing depth-confounded cells".
# build a null by repeatedly drawing a random (depth-agnostic) subsample of the FULL data at the same group sizes as the depth-matched analysis, and see how much disagreement with the full DGE that alone produces.
n_replicates <- 20

df_matched_size <- imap(list_sobj_matched,function(sobj,cell_id){
  tibble(cell_id = cell_id,
         n_pos = sum(sobj$TSPO_cat == "pos"),
         n_neg = sum(sobj$TSPO_cat == "neg"))
}) %>%
  bind_rows()

df_matched_size

random_subsample_DGE <- function(sobj,cell_id,n_pos,n_neg,seed){
  set.seed(seed)
  meta <- sobj@meta.data %>%
    rownames_to_column("barcodes")
  
  bc_pos <- meta %>% filter(TSPO_cat == "pos") %>% pull(barcodes) %>% sample(n_pos)
  bc_neg <- meta %>% filter(TSPO_cat == "neg") %>% pull(barcodes) %>% sample(n_neg)
  
  sobj_sub <- subset(sobj,cells = c(bc_pos,bc_neg))
  Idents(sobj_sub) <- "TSPO_cat"
  
  RunPresto(object = sobj_sub,ident.1 = "pos",ident.2 = "neg",min.pct = 0.01) %>%
    rownames_to_column("gene") %>%
    mutate(cell_id = cell_id)
}

df_DGE_random_null <- map(names(list_sobj),function(id){
  sizes <- df_matched_size %>% filter(cell_id == id)
  
  map(seq_len(n_replicates),function(rep_id){
    print(paste(id,"random null replicate",rep_id))
    random_subsample_DGE(list_sobj[[id]],id,
                         n_pos = sizes$n_pos,n_neg = sizes$n_neg,
                         seed = rep_id) %>%
      mutate(replicate = rep_id)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

df_DGE_random_null %>%
  write_tsv("../../out/table/custom/07_TSPO_res_test_posVSneg_randomSameSizeNull.tsv")

# compare each random replicate against the full DGE, using the same log2FC-correlation + DE-category-replication metrics used for the depth-matched comparison above
df_full_sig <- df_DGE_TSPO_full %>%
  mutate(DE_cat_full = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                                 avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                                 T~"no")) %>%
  dplyr::select(cell_id,gene,avg_log2FC_full = avg_log2FC,p_val_adj_full = p_val_adj,DE_cat_full)

df_random_compare <- df_DGE_random_null %>%
  dplyr::select(cell_id,gene,replicate,avg_log2FC_random = avg_log2FC,p_val_adj_random = p_val_adj) %>%
  mutate(DE_cat_random = case_when(avg_log2FC_random > 0.5 & p_val_adj_random < 0.01~"up",
                                   avg_log2FC_random < (-0.5) & p_val_adj_random < 0.01~"down",
                                   T~"no")) %>%
  inner_join(df_full_sig,by = c("cell_id","gene"))

df_random_compare %>%
  write_tsv("../../out/table/custom/07_TSPO_DGE_compare_full_vs_randomNull.tsv")

# per replicate: correlation of the effect size with the full analysis, and the fraction of full DEGs that stay significant + same-sign (replication rate)
df_random_null_summary <- df_random_compare %>%
  group_by(cell_id,replicate) %>%
  summarise(cor_log2FC_all = cor(avg_log2FC_full,avg_log2FC_random,method = "spearman"),
            .groups = "drop") %>%
  left_join(
    df_random_compare %>%
      filter(DE_cat_full != "no") %>%
      group_by(cell_id,replicate) %>%
      summarise(n_full_DE = n(),
                n_replicated = sum(DE_cat_random == DE_cat_full),
                replication_rate = n_replicated/n_full_DE,
                .groups = "drop"),
    by = c("cell_id","replicate")
  )

df_random_null_summary %>%
  write_tsv("../../out/table/custom/07_TSPO_DGE_randomNull_summary.tsv")

# the same two metrics for the depth-matched comparison, as a single observed value to place against the random same-size null distribution
df_matched_replication <- df_DGE_compare %>%
  group_by(cell_id) %>%
  summarise(cor_log2FC_all = cor(avg_log2FC_full,avg_log2FC_matched,method = "spearman"),
            .groups = "drop") %>%
  left_join(
    df_DGE_compare %>%
      filter(DE_cat_full != "no") %>%
      group_by(cell_id) %>%
      summarise(n_full_DE = n(),
                n_replicated = sum(DE_cat_matched == DE_cat_full),
                replication_rate = n_replicated/n_full_DE,
                .groups = "drop"),
    by = "cell_id"
  )

df_matched_replication %>%
  write_tsv("../../out/table/custom/07_TSPO_DGE_depthMatched_replicationRate.tsv")

# plot: black = random same-size (depth-agnostic) null distribution, red diamond = the actual depth-matched result. if the red diamond falls within the null cloud, the drop from the full result is explained by the smaller n alone; if it falls clearly below the null, depth-associated cells were carrying real signal
p_repl <- df_random_null_summary %>%
  ggplot(aes(x = cell_id,y = replication_rate)) +
  geom_jitter(width = 0.1,alpha = 0.5,size = 1) +
  geom_point(data = df_matched_replication,aes(x = cell_id,y = replication_rate),
             col = "red",size = 3,shape = 18) +
  ylim(0,1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(y = "replication rate of full DEGs (same sign + significant)")

p_cor <- df_random_null_summary %>%
  ggplot(aes(x = cell_id,y = cor_log2FC_all)) +
  geom_jitter(width = 0.1,alpha = 0.5,size = 1) +
  geom_point(data = df_matched_replication,aes(x = cell_id,y = cor_log2FC_all),
             col = "red",size = 3,shape = 18) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(y = "avg_log2FC spearman correlation with full")

p_repl + p_cor
ggsave("../../out/image/custom/07_TSPO_randomNull_vs_depthMatched.pdf",width = 10,height = 5)

# pathway-level comparison: full vs depth-matched vs random-null DGE --------
# gene-level concordance can be noisy for individual genes sitting near the significance boundary; aggregating into pathways is a more stable way to ask whether the *biological interpretation* survives depth matching. doing it against the same random same-size null used above separates "pathway signal lost from smaller n" from "pathway signal lost from removing depth-confounded cells"
library(fgsea)
library(msigdbr)

gene_sets_hallmark <- msigdbr(species = "Homo sapiens",category = "H")
pathways_hallmark <- split(x = gene_sets_hallmark$gene_symbol,f = gene_sets_hallmark$gs_name)

build_ranks <- function(df){
  df_rank <- df %>%
    filter(!is.na(avg_log2FC)) %>%
    group_by(gene) %>%
    summarise(avg_log2FC = mean(avg_log2FC),.groups = "drop")

  setNames(df_rank$avg_log2FC,df_rank$gene)
}

run_GSEA_TSPO <- function(ranks,seed = 123){
  set.seed(seed)
  fgsea(pathways_hallmark,ranks,minSize = 10,maxSize = 500) %>%
    as_tibble()
}

# GSEA on the full and depth-matched DGE, per cell type
df_GSEA_full <- map(names(list_sobj),function(id){
  ranks <- df_DGE_TSPO_full %>% filter(cell_id == id) %>% build_ranks()
  run_GSEA_TSPO(ranks) %>% mutate(cell_id = id)
}) %>%
  bind_rows()

df_GSEA_matched <- map(names(list_sobj),function(id){
  ranks <- df_DGE_TSPO_matched %>% filter(cell_id == id) %>% build_ranks()
  run_GSEA_TSPO(ranks) %>% mutate(cell_id = id)
}) %>%
  bind_rows()

# GSEA on each random same-size null replicate, per cell type
df_GSEA_random_null <- map(names(list_sobj),function(id){
  map(seq_len(n_replicates),function(rep_id){
    ranks <- df_DGE_random_null %>%
      filter(cell_id == id,replicate == rep_id) %>%
      build_ranks()
    run_GSEA_TSPO(ranks,seed = rep_id) %>%
      mutate(cell_id = id,replicate = rep_id)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

# save the full pathway-level tables (leadingEdge collapsed to a string to save as tsv)
bind_rows(df_GSEA_full %>% mutate(condition = "full"),
          df_GSEA_matched %>% mutate(condition = "matched")) %>%
  mutate(leadingEdge = unlist(lapply(.$leadingEdge,function(x){paste0(x,collapse = "|")}))) %>%
  arrange(cell_id,condition,padj,-abs(NES)) %>%
  write_tsv("../../out/table/custom/07_TSPO_GSEA_hallmark_full_vs_depthMatched.tsv")

df_GSEA_random_null %>%
  mutate(leadingEdge = unlist(lapply(.$leadingEdge,function(x){paste0(x,collapse = "|")}))) %>%
  write_tsv("../../out/table/custom/07_TSPO_GSEA_hallmark_randomSameSizeNull.tsv")

# full vs matched: NES correlation and replication rate (mirrors the gene-level check)
df_GSEA_compare_matched <- df_GSEA_full %>%
  dplyr::select(cell_id,pathway,NES_full = NES,padj_full = padj) %>%
  inner_join(df_GSEA_matched %>%
               dplyr::select(cell_id,pathway,NES_matched = NES,padj_matched = padj),
             by = c("cell_id","pathway"))

df_GSEA_matched_replication <- df_GSEA_compare_matched %>%
  group_by(cell_id) %>%
  summarise(cor_NES = cor(NES_full,NES_matched,method = "spearman"),
            .groups = "drop") %>%
  left_join(
    df_GSEA_compare_matched %>%
      filter(padj_full < 0.05) %>%
      group_by(cell_id) %>%
      summarise(n_full_sig = n(),
                n_replicated = sum(padj_matched < 0.05 & sign(NES_full) == sign(NES_matched)),
                replication_rate = n_replicated/n_full_sig,
                .groups = "drop"),
    by = "cell_id"
  )

df_GSEA_matched_replication %>%
  write_tsv("../../out/table/custom/07_TSPO_GSEA_depthMatched_replicationRate.tsv")

# full vs each random same-size replicate: same two metrics, building the null distribution
df_GSEA_compare_random <- df_GSEA_full %>%
  dplyr::select(cell_id,pathway,NES_full = NES,padj_full = padj) %>%
  inner_join(df_GSEA_random_null %>%
               dplyr::select(cell_id,pathway,replicate,NES_random = NES,padj_random = padj),
             by = c("cell_id","pathway"))

df_GSEA_random_null_summary <- df_GSEA_compare_random %>%
  group_by(cell_id,replicate) %>%
  summarise(cor_NES = cor(NES_full,NES_random,method = "spearman"),
            .groups = "drop") %>%
  left_join(
    df_GSEA_compare_random %>%
      filter(padj_full < 0.05) %>%
      group_by(cell_id,replicate) %>%
      summarise(n_full_sig = n(),
                n_replicated = sum(padj_random < 0.05 & sign(NES_full) == sign(NES_random)),
                replication_rate = n_replicated/n_full_sig,
                .groups = "drop"),
    by = c("cell_id","replicate")
  )

df_GSEA_random_null_summary %>%
  write_tsv("../../out/table/custom/07_TSPO_GSEA_randomNull_summary.tsv")

# plot: same logic as the gene-level check, black = random same-size null,
# red diamond = depth-matched, now at the pathway level
p_GSEA_repl <- df_GSEA_random_null_summary %>%
  ggplot(aes(x = cell_id,y = replication_rate)) +
  geom_jitter(width = 0.1,alpha = 0.5,size = 1) +
  geom_point(data = df_GSEA_matched_replication,aes(x = cell_id,y = replication_rate),
             col = "red",size = 3,shape = 18) +
  ylim(0,1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(y = "replication rate of full sig. Hallmark pathways (same sign + significant)")

p_GSEA_cor <- df_GSEA_random_null_summary %>%
  ggplot(aes(x = cell_id,y = cor_NES)) +
  geom_jitter(width = 0.1,alpha = 0.5,size = 1) +
  geom_point(data = df_GSEA_matched_replication,aes(x = cell_id,y = cor_NES),
             col = "red",size = 3,shape = 18) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(y = "NES spearman correlation with full")

p_GSEA_repl + p_GSEA_cor
ggsave("../../out/image/custom/07_TSPO_GSEA_randomNull_vs_depthMatched.pdf",width = 10,height = 5)

# scatter of NES for pathways significant in either full or depth-matched
df_GSEA_compare_matched %>%
  mutate(sig_either = padj_full < 0.05 | padj_matched < 0.05) %>%
  filter(sig_either) %>%
  ggplot(aes(x = NES_full,y = NES_matched)) +
  geom_point(alpha = 0.5,size = 1) +
  ggrepel::geom_text_repel(aes(label = pathway),size = 1.8,max.overlaps = 15) +
  geom_abline(slope = 1,intercept = 0,linetype = "dashed",col = "red") +
  facet_wrap(~cell_id,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/custom/07_TSPO_GSEA_scatter_NES_full_vs_depthMatched.pdf",width = 10,height = 5)

# volcano comparison: full vs average random-null simulation vs depth-matched ----
# summarise the 20 random same-size replicates into one representative "average simulation" estimate per gene (mean effect size across replicates; median p_val_adj, more robust to a single noisy replicate than the mean), so it can be plotted side by side with the full and the depth-matched estimates on the same volcano layout
df_DGE_random_avg <- df_DGE_random_null %>%
  group_by(cell_id,gene) %>%
  summarise(avg_log2FC = mean(avg_log2FC),
            p_val_adj = median(p_val_adj),
            .groups = "drop")

df_volcano_compare <- bind_rows(
  df_DGE_TSPO_full %>% dplyr::select(cell_id,gene,avg_log2FC,p_val_adj) %>% mutate(condition = "full"),
  df_DGE_random_avg %>% mutate(condition = "avg_simulation"),
  df_DGE_TSPO_matched %>% dplyr::select(cell_id,gene,avg_log2FC,p_val_adj) %>% mutate(condition = "depth_matched")
) %>%
  mutate(condition = factor(condition,levels = c("full","avg_simulation","depth_matched"))) %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no"))

df_volcano_compare %>%
  write_tsv("../../out/table/custom/07_TSPO_volcano_compare_full_avgSimulation_depthMatched.tsv")

df_volcano_compare %>%
  filter(gene == "TREM2")

ggplot() +
  geom_point(data = df_volcano_compare %>% filter(DE_cat == "no"),
             aes(x = avg_log2FC,y = -log10(p_val_adj)),size = 0.5,alpha = 0.2) +
  geom_point(data = df_volcano_compare %>% filter(DE_cat != "no"),
             aes(x = avg_log2FC,y = -log10(p_val_adj)),size = 0.5,alpha = 0.4,col = "red") +
  ggrepel::geom_text_repel(data = df_volcano_compare %>%
                              filter(DE_cat != "no") %>%
                              group_by(cell_id,condition) %>%
                              slice_min(order_by = p_val_adj,n = 15),
                            aes(x = avg_log2FC,y = -log10(p_val_adj),label = gene),
                            size = 2,max.overlaps = 15) +
  facet_grid(cell_id~condition,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/custom/07_TSPO_volcano_full_vs_avgSimulation_vs_depthMatched.pdf",width = 15,height = 10)

# enrichR on the DE gene lists: full vs avg_simulation vs depth_matched -----
# same DB panel/workflow as in 40_plotGene_integrationSkip_manualClean_harmony_update.R, run separately on the up/down gene lists of every cell type x condition, to check whether the pathway-level interpretation still holds regardless of the QC-driven concern
library(enrichR)
library(scales)

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_Pathways_2024","HDSigDB_Human_2021","Azimuth_Cell_Types_2021","GO_Biological_Process_2023","Descartes_Cell_Types_and_Tissue_2021","CellMarker_Augmented_2021")

# pull the up/down gene lists for every cell type x condition combination
# (a group_id only exists if it has at least one DE gene, so empty lists are
# skipped automatically)
list_genes_TSPO <- df_volcano_compare %>%
  filter(DE_cat != "no") %>%
  mutate(group_id = paste(cell_id,condition,DE_cat,sep = "|")) %>%
  {split(.$gene,.$group_id)}

list_enrichr_TSPO <- imap(list_genes_TSPO,function(genes,group_id){
  print(group_id)
  out_enrich <- enrichr(genes,dbs_db)

  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>%
    unlist()

  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation") %>%
    mutate(group_id = group_id)
}) %>%
  bind_rows() %>%
  separate(group_id,into = c("cell_id","condition","DE_cat"),sep = "\\|")

list_enrichr_TSPO %>%
  write_tsv("../../out/table/custom/07_TSPO_enrichR_full_vs_avgSimulation_vs_depthMatched.tsv")

# plot top 10 terms per db, per cell type x condition x direction
plot_list_TSPO <- list_enrichr_TSPO %>%
  mutate(facet_id = paste(cell_id,condition,DE_cat,sep = " | ")) %>%
  split(f = .$facet_id)

list_plot_TSPO <- imap(plot_list_TSPO,function(x,facet_id){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term,Combined.Score,.desc = F)) %>%
    ggplot(aes(y = Term,x = Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) +
    geom_point() +
    facet_wrap(~annotation,scales = "free",ncol = 1) +
    theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2)) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    ggtitle(facet_id)
})

# one panel per cell type, combining full / avg_simulation / depth_matched x up/down
walk(names(list_sobj),function(id){
  wrap_plots(list_plot_TSPO[str_starts(names(list_plot_TSPO),id)])
  ggsave(paste0("../../out/image/custom/07_TSPO_enrichR_",id,"_full_vs_avgSimulation_vs_depthMatched.pdf"),
         width = 20,height = 25,limitsize = FALSE)
})

# correlation of the enrichR score across the three approaches --------------
# use Combined.Score (magnitude x significance of the overrepresentation) as the enrichment score, and compare it term-by-term (db x pathway x direction) between full, avg_simulation and depth_matched
df_enrichr_wide <- list_enrichr_TSPO %>%
  dplyr::select(cell_id,condition,DE_cat,annotation,Term,Combined.Score,Adjusted.P.value) %>%
  pivot_wider(names_from = condition,values_from = c(Combined.Score,Adjusted.P.value))

df_enrichr_wide %>%
  write_tsv("../../out/table/custom/07_TSPO_enrichR_CombinedScore_wide.tsv")

df_enrichr_cor <- df_enrichr_wide %>%
  group_by(cell_id) %>%
  summarise(cor_full_vs_matched = cor(Combined.Score_full,Combined.Score_depth_matched,method = "spearman",use = "pairwise.complete.obs"),
            cor_full_vs_avgSim = cor(Combined.Score_full,Combined.Score_avg_simulation,method = "spearman",use = "pairwise.complete.obs"),
            cor_avgSim_vs_matched = cor(Combined.Score_avg_simulation,Combined.Score_depth_matched,method = "spearman",use = "pairwise.complete.obs"),
            .groups = "drop")

df_enrichr_cor %>%
  write_tsv("../../out/table/custom/07_TSPO_enrichR_CombinedScore_correlation.tsv")

df_enrichr_cor

# scatter of Combined.Score for the three pairwise comparisons, restricted
# to terms significant (Adjusted.P.value<0.05) in at least one of the two
# approaches being compared, for readability
df_enrichr_pairs <- bind_rows(
  df_enrichr_wide %>%
    transmute(cell_id,DE_cat,annotation,Term,
              score_x = Combined.Score_full,score_y = Combined.Score_depth_matched,
              sig = Adjusted.P.value_full < 0.05 | Adjusted.P.value_depth_matched < 0.05,
              comparison = "full vs depth_matched"),
  df_enrichr_wide %>%
    transmute(cell_id,DE_cat,annotation,Term,
              score_x = Combined.Score_full,score_y = Combined.Score_avg_simulation,
              sig = Adjusted.P.value_full < 0.05 | Adjusted.P.value_avg_simulation < 0.05,
              comparison = "full vs avg_simulation"),
  df_enrichr_wide %>%
    transmute(cell_id,DE_cat,annotation,Term,
              score_x = Combined.Score_avg_simulation,score_y = Combined.Score_depth_matched,
              sig = Adjusted.P.value_avg_simulation < 0.05 | Adjusted.P.value_depth_matched < 0.05,
              comparison = "avg_simulation vs depth_matched")
) %>%
  mutate(comparison = factor(comparison,levels = c("full vs avg_simulation","full vs depth_matched","avg_simulation vs depth_matched")))

df_enrichr_pairs %>%
  filter(sig) %>%
  ggplot(aes(x = score_x,y = score_y)) +
  geom_point(alpha = 0.4,size = 0.8) +
  geom_abline(slope = 1,intercept = 0,linetype = "dashed",col = "red") +
  facet_grid(cell_id~comparison,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  labs(x = "Combined.Score_full",y = "Combined.Score_test")
ggsave("../../out/image/custom/07_TSPO_enrichR_scatter_CombinedScore_compare.pdf",width = 14,height = 8)

# TSPO DGE pos vs neg split by disease (CTRL/MS) ---------------------------
# same DGE as in 40_..., stratified within CTRL and within MS
run_DGE_TSPO_disease <- function(sobj,cell_id){
  map(sort(unique(sobj$disease)),function(disease_id){
    print(paste(cell_id,disease_id))
    sobj_sub <- subset(sobj,subset = disease == disease_id)
    Idents(sobj_sub) <- "TSPO_cat"
    
    table(Idents(sobj_sub))
    
    # avg_log2FC: positive values indicate higher expression in the pos group
    RunPresto(object = sobj_sub,ident.1 = "pos",ident.2 = "neg") %>%
      rownames_to_column("gene") %>%
      mutate(cell_id = cell_id,disease = disease_id)
  }) %>%
    bind_rows()
}

df_DGE_TSPO_disease <- imap(list_sobj,run_DGE_TSPO_disease) %>%
  bind_rows()

df_DGE_TSPO_disease %>%
  write_tsv("../../out/table/custom/07_TSPO_res_test_posVSneg_byDisease.tsv")

# annotate DE genes and plot a volcano per cell type / disease
df_volcano_TSPO_disease <- df_DGE_TSPO_disease %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no"))

ggplot() +
  geom_point(data = df_volcano_TSPO_disease %>% filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=0.5,alpha=0.2) +
  geom_point(data = df_volcano_TSPO_disease %>% filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=0.5,alpha=0.2,col="red") +
  ggrepel::geom_text_repel(data = df_volcano_TSPO_disease %>% filter(DE_cat!="no") %>% group_by(cell_id,disease) %>% slice_min(order_by = p_val_adj,n = 15),
                           aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene),size = 2,max.overlaps = 15) +
  facet_wrap(cell_id~disease,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/custom/07_TSPO_volcano_test_posVSneg_byDisease.pdf",width = 12,height = 10)
