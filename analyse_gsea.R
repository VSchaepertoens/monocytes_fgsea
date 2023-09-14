## 
##
## Script name: monocytes_gsea
##
## Purpose of script: runs fGSEA analysis on three datasets (alwof_vs_hpyl, hpyl_vs_lps, alwof_vs_lps)
##
## Author: Dr. Veronika Schäpertöns & Prof. Nikolaus Fortelny
##
## Date Created: 23.03.2023
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


library(readxl)
library(tidyverse)
library(data.table)
library(janitor)
library(fgsea)
library(fs)



# import data -------------------------------------------------------------

coef_names <- c(
  t1lp_svs_control = "lps_vs_uninduced",
  t2hpylvs_control = "hpyl_vs_uninduced",
  t3alwofvs_control = "alwof_vs_uninduced",
  t2hpylvs_t1lps = "hpyl_vs_lps",
  t3alwofvs_t1lps = "alwof_vs_lps",
  t3alwofvs_t2hpyl = "alwof_vs_hpyl"
)

identifiers <-
  read_csv("metadata/identifiers.csv") %>%
  select(accession, gene_names)

df_data <- 
  read_csv("analysis/results_limma.csv") %>%
  clean_names() %>% 
  select(accession, starts_with("coef"), starts_with("p_value_adj")) %>% 
  pivot_longer(
    !accession,
    names_to = c("type", "comparison"),
    names_pattern = "(p_value_adj|coef)_(.+)"
  ) %>%
  pivot_wider(names_from = type) %>%
  rename(log_fc = coef, padj = p_value_adj, coef = comparison) %>%
  left_join(identifiers, by = "accession") %>% 
  mutate(
    coef = recode(coef, !!!coef_names),
    gene = str_replace(gene_names, " .+$", "")
  )

df_data


# load genesets (pathways for fgsea) ------------------------------

genesets <- read_tsv("analysis/enrichr_database.tsv")


# enrichment analysis using fgsea -----------------------------------------
set.seed(119)
gsea.res <- data.table()
de.grp <- df_data$coef[1]

for (de.grp in unique(df_data$coef)) {
  print(de.grp)
  gsea.res <- rbind(gsea.res, 
                    data.table(fgsea(
                      pathways = with(genesets, split(Gene, paste(DB, Geneset, sep = "__"))),
                      stats = with(df_data[df_data$coef == de.grp,], setNames(object = log_fc, nm = gene)),
                      nperm = 1e6), 
                    grp = de.grp))
}

gsea.res_orig <- gsea.res
#gsea.res <- gsea.res_orig

gsea.res$leadingEdge <- sapply(
  gsea.res$leadingEdge, 
  function(vec) paste(vec, collapse = ",")
)


# gsea results ------------------------------------------------------------

dim(gsea.res[padj < 0.05])

gsea.res[padj < 0.05]
gsea.res[padj < 0.05][,-"leadingEdge"][order(padj)]

gsea.res[padj < 0.05][,-"leadingEdge"][grp == "hpyl_vs_lps"][order(padj)]
gsea.res[padj < 0.05][,-"leadingEdge"][grp == "alwof_vs_lps"][order(padj)]
gsea.res[padj < 0.05][,-"leadingEdge"][grp == "alwof_vs_hpyl"][order(padj)]
gsea.res[padj < 0.05][,-"leadingEdge"][grp == "alwof_vs_uninduced"][order(padj)]
gsea.res[padj < 0.05][,-"leadingEdge"][grp == "hpyl_vs_uninduced"][order(padj)]
gsea.res[padj < 0.05][,-"leadingEdge"][grp == "lps_vs_uninduced"][order(padj)]


dir_create("analysis")
write_tsv(x = gsea.res, file = "analysis/results_gsea.tsv")

