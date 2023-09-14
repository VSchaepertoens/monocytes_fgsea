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


# make plots --------------------------------------------------------------

plot_gsea <- function(grp, plot_title, revert = FALSE) {
  data_to_plot <- 
    gsea.res %>% 
    filter(padj < 0.05, grp == {{grp}}) %>% 
    arrange(padj)
  
  if (revert)
    data_to_plot <- 
      data_to_plot %>% 
      mutate(NES = -NES)
  
  ggplot(data_to_plot) +
    geom_col(
      mapping = aes(x = NES,
                    y = reorder(pathway_short, NES),
                    fill = -log10(padj))
    ) +
    xlab("normalized enrichment score") +
    ylab("gene sets") +
    ggtitle(plot_title) +
    scale_fill_gradient(
      low = "black",
      high = "red",
      name = expression(-log[10](p[adj]))
    ) +
    theme_bw() +
    theme(text = element_text(size = 10))
}

# aci vs hp
plot_gsea("alwof_vs_hpyl", "Enriched in Hp vs Aci", revert = TRUE)

ggsave("figures/gsea_hp_vs_aci_600dpi.png",
       height = 150,
       width = 300,
       units = "mm",
       dpi = 600)


# lps vs hp
plot_gsea("hpyl_vs_lps", "Enriched in Hp vs LPS")

ggsave("figures/gsea_hp_vs_lps_600dpi.png",
       height = 200,
       width = 300,
       units = "mm",
       dpi = 600)