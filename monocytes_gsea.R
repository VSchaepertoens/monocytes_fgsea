## ---------------------------
##
## Script name: monocytes_gsea
##
## Purpose of script: runs fGSEA analysis on three datasets (alwof_vs_hpyl, hpyl_vs_lps, alwof_vs_lps)
##
## Author: Dr. Veronika Schäpertöns & Nikolaus Fortelny
##
## Date Created: 23.03.2023
##
## Copyright (c) Veronika Schäpertöns, 2023
## Email: veronika.schaepertoens@plus.ac.at
##
## ---------------------------
##
## Notes:
##   
##
## 


# import libraries ---------------------------------------------------------
library(readxl)
library(tidyverse)
library(data.table)
library(janitor)
library(fgsea)

# import data -------------------------------------------------------------

df_data_hpyl_vs_lps <- read_xlsx("data/Monocytes_Experiment_Results_20230316.xlsx", sheet = "Hpyl_vs_LPS") %>%
  clean_names() %>% 
  rename(log_fc = coef_t2hpylvs_t1lps, 
         padj = p_value_adj_t2hpylvs_t1lps) %>%
  mutate(coef = "hpyl_vs_lps")

df_data_alwof_vs_lps <- read_xlsx("data/Monocytes_Experiment_Results_20230316.xlsx", sheet = "Alwof_vs_LPS") %>%
  clean_names() %>% 
  rename(log_fc = coef_t3alwofvs_t1lps, 
         padj = p_value_adj_t3alwofvs_t1lps) %>%
  mutate(coef = "alwof_vs_lps")

df_data_alwof_vs_hpyl <- read_xlsx("data/Monocytes_Experiment_Results_20230316.xlsx", sheet = "AlwofvsHpyl") %>%
  clean_names() %>% 
  rename(log_fc = coef_t3alwofvs_t2hpyl, 
         padj = p_value_adj_t3alwofvs_t2hpyl) %>%
  mutate(coef = "alwof_vs_hpyl")  


# bind dataframes by rows & clean gene names -------------------------------------------------------

df_data <- bind_rows(
  df_data_alwof_vs_hpyl,
  df_data_alwof_vs_lps,
  df_data_hpyl_vs_lps
)

df_data['gene'] <- gsub(" .+$", "", df_data$gene_names)


# load genesets(actually pathways for fgsea) ------------------------------

genesets <- read_tsv("database/enrichr/enrichr_database.tsv")


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

gsea.res$leadingEdge <- sapply(gsea.res$leadingEdge, function(vec) paste(vec, collapse = ","))

dim(gsea.res[padj < 0.05])

gsea.res[padj < 0.05]

gsea.res[padj < 0.05][,-"leadingEdge"][order(padj)]

gsea.res[padj < 0.05][,-"leadingEdge"][grp == "hpyl_vs_lps"][order(padj)]
gsea.res[padj < 0.05][,-"leadingEdge"][grp == "alwof_vs_lps"][order(padj)]
gsea.res[padj < 0.05][,-"leadingEdge"][grp == "alwof_vs_hpyl"][order(padj)]


write_tsv(x = gsea.res, file = "analysis/gsea_alwof_vs_hpyl.tsv")

