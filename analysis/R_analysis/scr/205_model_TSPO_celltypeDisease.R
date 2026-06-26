# AIM ---------------------------------------------------------------------
# model TSPO by ell type and disease

library(readr)
library(dplyr)
library(ggplot2)
library(emmeans)
library(scales)

df_avg_cellType <- read_tsv("../../out/table/revision/df_avg_cellType_TSPO_update.tsv") %>%
  filter(!(sample %in% c("s4","s26","s22")))

head(df_avg_cellType)
table(df_avg_cellType$expertAnno.l1)
nrow(df_avg_cellType)

test_lm_cellType_disease <- lm(avg_exp ~ expertAnno.l1 * disease, data = df_avg_cellType)
summary(test_lm_cellType_disease)

emm_test_lm_cellType_disease <- emmeans(test_lm_cellType_disease, ~ disease | expertAnno.l1)
pairs(emm_test_lm_cellType_disease, reverse = T)

df_emm2 <- data.frame(emm_test_lm_cellType_disease)

df_emm2$expertAnno.l1 <- factor(
  df_emm2$expertAnno.l1,
  levels = c("AST", "EPENDYMA" ,"LYM", "IMM", "EXC.NEU", "INH.NEU", "OLIGO", "OPC", "VAS"))

head(df_emm2)

ggplot() +
  geom_point(data = df_avg_cellType, aes(x = expertAnno.l1, y = avg_exp,col=disease),
             alpha = 0.5, shape = 1,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0,dodge.width = 0.8)) +
  geom_point(data = df_emm2, aes(x = expertAnno.l1, y = emmean,col=disease),
             position = position_dodge(width = 0.8),
             size = 3.5, shape = 18) +
  geom_errorbar(data = df_emm2,
                aes(x = expertAnno.l1, ymin = lower.CL, ymax = upper.CL,col=disease),
                position = position_dodge(width = 0.8),
                width = 0.15, linewidth = 0.7) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    # legend.position = "none",
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Cell type",
    y = "TSPO average expression"
  ) +
  coord_cartesian(ylim = c(0, 0.4)) +
  scale_x_discrete(labels = c(
    "ASTRO" = "Astrocytes",
    "LYM" = "Lymphocytes",
    "MG" = "Microglia",
    "NEU_EXC" = "Excitatory neurons",
    "NEU_INH" = "Inhibitory neurons",
    "NEU_PY" = "Pyramidal neurons",
    "OLIG" = "Oligodendrocytes",
    "OPC" = "OPC",
    "VAS" = "Vascular cells"
  )) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 5)) +
  ggtitle("TSPO expression by cell type")


ggplot() +
  geom_boxplot(data = df_avg_cellType, aes(x = expertAnno.l1, y = avg_exp,col=disease),outlier.shape = NA) +
  geom_point(data = df_avg_cellType, aes(x = expertAnno.l1, y = avg_exp,col=disease),
             alpha = 0.5, shape = 1,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0,dodge.width = 0.8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    # legend.position = "none",
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 0.4)) +
  scale_x_discrete(labels = c(
    "ASTRO" = "Astrocytes",
    "LYM" = "Lymphocytes",
    "MG" = "Microglia",
    "NEU_EXC" = "Excitatory neurons",
    "NEU_INH" = "Inhibitory neurons",
    "NEU_PY" = "Pyramidal neurons",
    "OLIG" = "Oligodendrocytes",
    "OPC" = "OPC",
    "VAS" = "Vascular cells"
  )) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 5)) +
  ggtitle("TSPO expression by cell type")
