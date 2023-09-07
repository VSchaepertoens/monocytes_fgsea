## ---------------------------
##
## Script name: plot_gsea
##
## Purpose of script: plots NES of results from the fGSEA analysis
##
## Author: Dr. Veronika Schäpertöns
##
## Date Created: 07.09.2023
##
## Copyright (c) Veronika Schäpertöns, 2023
## Email: veronika.schaepertoens@plus.ac.at
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(ggplot2)
library(tidyverse)
library(data.table)
library(paletteer)


gsea.res <- read_tsv("analysis/results_gsea.tsv") %>%
  data.table()

dim(gsea.res[padj < 0.05])

gsea.res[padj < 0.05]

# aci vs hp ---------------------------------------------------------------

data_to_plot <- gsea.res[padj < 0.05][,-"leadingEdge"][grp == "alwof_vs_hpyl"][order(padj)]

write.table(data_to_plot, "analysis/gsea_alwof_hpyl.tsv")

data_to_plot <- read_delim("analysis/gsea_alwof_hpyl_short.tsv",
                           delim = " ")
## bar graph ---------------------------------------------------------------

png("figures/gsea_hp_vs_aci_600dpi.png",
    height = 150,
    width = 150,
    unit = "mm",
    res = 600)

ggplot(data = data_to_plot) +
  geom_col(
    mapping = aes(x = (NES*-1),
                  y = reorder(pathway, NES),
                  fill = -log10(padj))
  ) +
  scale_fill_paletteer_c("ggthemes::Classic Red-Black", direction = -1) +
  xlab("normalized enrichment score") +
  ylab("gene sets") +
  ggtitle("Enriched in Hp vs Aci") +
  scale_fill_gradient(low = "black", high = "red", name = expression(-log[10](p[adj]))) +
  theme_bw() +
  theme(text = element_text(size = 10, family = "sans"))

dev.off()


# lps vs hp ---------------------------------------------------------------

data_to_plot <- gsea.res[padj < 0.05][,-"leadingEdge"][grp == "hpyl_vs_lps"][order(padj)]

write.table(data_to_plot, "analysis/gsea_hpyl_lps.tsv")

data_to_plot <- read_delim("analysis/gsea_hpyl_lps_short.tsv",
                           delim = " ")
## bar graph ---------------------------------------------------------------

png("figures/gsea_hp_vs_lps_600dpi.png",
    height = 200,
    width = 200,
    unit = "mm",
    res = 600)

ggplot(data = data_to_plot) +
  geom_col(
    mapping = aes(x = (NES),
                  y = reorder(pathway, NES),
                  fill = -log10(padj))
  ) +
  scale_fill_paletteer_c("ggthemes::Classic Red-Black", direction = -1) +
  xlab("normalized enrichment score") +
  ylab("gene sets") +
  ggtitle("Enriched in Hp vs LPS") +
  scale_fill_gradient(low = "black", high = "red", name = expression(-log[10](p[adj]))) +
  theme_bw() +
  theme(text = element_text(size = 10, family = "sans"))

dev.off()
