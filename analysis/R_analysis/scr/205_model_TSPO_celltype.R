# AIM ---------------------------------------------------------------------
# test model TSPO per cell type

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
test_lm_cellType_lin <- lm(avg_exp ~ expertAnno.l1, data = df_avg_cellType)
summary(test_lm_cellType_lin)

emm_test_lm_cellType_lin <- emmeans(test_lm_cellType_lin, ~ expertAnno.l1)
pairs(emm_test_lm_cellType_lin, reverse = T)

df_emm2 <- as.data.frame(emm_test_lm_cellType_lin)

df_emm2$expertAnno.l1 <- factor(
  df_emm2$expertAnno.l1,
  levels = c("AST", "EPENDYMA" ,"LYM", "IMM", "EXC.NEU", "INH.NEU", "OLIGO", "OPC", "VAS"))

head(df_emm2)

ggplot() +
  geom_point(data = df_avg_cellType, aes(x = expertAnno.l1, y = avg_exp),
             alpha = 0.8, shape = 1,
             position = position_jitter(width = 0.2, height = 0)) +
  geom_point(data = df_emm2, aes(x = expertAnno.l1, y = emmean),
             size = 3, shape = 18) +
  geom_errorbar(data = df_emm2,
                aes(x = expertAnno.l1, ymin = lower.CL, ymax = upper.CL),
                width = 0.15, linewidth = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "none",
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Cell type",
    y = "TSPO average expression"
  ) +
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
  geom_point(data = df_avg_cellType, aes(x = expertAnno.l1, y = avg_exp),
             alpha = 0.5, shape = 1, color = "grey60",
             position = position_jitter(width = 0.2, height = 0)) +
  geom_point(data = df_emm2, aes(x = expertAnno.l1, y = emmean),
             size = 3.5, shape = 18, color = "red") +
  geom_errorbar(data = df_emm2,
                aes(x = expertAnno.l1, ymin = lower.CL, ymax = upper.CL),
                width = 0.15, linewidth = 0.7, color = "red") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "none",
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
