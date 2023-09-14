## 
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
## 
##
## Notes:
##   
##
## 

library(tidyverse)



# load data ---------------------------------------------------------------

database_short_names <- c(
  "CORUM" = "CORUM",
  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X" = "ENCODE",
  "GO_Biological_Process_2018" = "GO",
  "GO_Cellular_Component_2018" = "GO",
  "GO_Molecular_Function_2018" = "GO",
  "Human_Gene_Atlas" = "HGA",
  "KEGG_2019_Human" = "KEGG",
  "MSigDB_Hallmark_2020" = "MSigDB",
  "NCI-Nature_2016" = "NCI",
  "TRANSFAC_and_JASPAR_PWMs" = "TRANSFAC",
  "TRRUST_Transcription_Factors_2019" = "TRRUST",
  "WikiPathways_2019_Human" = "WikiPathways_2019_Human"
)

pathway_short_names <- 
  read_tsv("analysis/results_gsea.tsv") %>% 
  distinct(pathway) %>% 
  separate_wider_delim(
    pathway,
    delim = "__",
    names = c("database", "pathway_short"),
    cols_remove = FALSE
  ) %>% 
  mutate(
    database_short = recode(database, !!!database_short_names),
    pathway_short = 
      pathway_short %>% 
      str_replace(" \\(GO:\\d+\\)", "") %>% 
      str_replace(" WP\\d+", "")
  ) %>% 
  unite(pathway_short, database_short, col = "pathway_short", sep = " ") %>% 
  select(pathway, pathway_short) %>% 
  mutate(pathway_short = make.unique(pathway_short))

gsea.res <-
  read_tsv("analysis/results_gsea.tsv") %>% 
  left_join(pathway_short_names, by = "pathway")



# bar graph ---------------------------------------------------------------

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
# bar graph ---------------------------------------------------------------

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
