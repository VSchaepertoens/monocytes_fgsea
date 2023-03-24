#GSEA analysis for monocytes data induced with hpylori/alwoffii/LPS


# import libraries ---------------------------------------------------------
library(readxl)
library(tidyverse)
library(janitor)
library(fgsea)

# import data -------------------------------------------------------------

df_data <- read_xlsx("data/Monocytes_Experiment_Results_20230316.xlsx", sheet = 5)%>%
  clean_names()


# clean gene names -------------------------------------------------------

df_data['gene'] <- gsub(" .+$", "", df_data$gene_names)


# load genesets(actually pathways for fgsea) ------------------------------

genesets <- read_tsv("database/enrichr/enrichr_database.tsv")


# enrichment analysis using fgsea -----------------------------------------
set.seed(119)

gsea.res <- data.table(fgsea(
  pathways = with(genesets, split(Gene, paste(DB, Geneset, sep="__"))),
  stats = with(df_data, setNames(object = df_data$coef_t3alwofvs_t2hpyl, nm = df_data$gene)),
  nperm=1e6
  )
)

gsea.res$leadingEdge <- sapply(gsea.res$leadingEdge, function(vec) paste(vec, collapse = ","))

dim(gsea.res[padj < 0.05])

gsea.res[padj < 0.05]

gsea.res[padj < 0.05][,-"leadingEdge"][order(padj)]

write_tsv(x = gsea.res, file = "analysis/gsea_alwof_vs_hpyl.tsv")

