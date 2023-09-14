library(tidyverse)
library(janitor)
library(pheatmap)
library(limma)
library(fs)

dir_create("figures")



# load data ---------------------------------------------------------------

data_matrix <-
  read_delim(
    "data/proteinGroups_log2_substractMedian_Batchcorrected.csv",
    delim = ";"
  ) %>% 
  column_to_rownames(var = "Accession") %>%
  as.matrix()

colnames(data_matrix) <- str_replace_all(
  colnames(data_matrix),
  c(Hpyl = "Hp", Alwof = "Aci", uninduced = "Uninduced")
)

df_data_unind_vs_treated <-
  read_csv("analysis/results_limma.csv") %>%
  clean_names()



# data exploration PCA & correlations -------------------------------------

plot(density(data_matrix))

boxplot(data_matrix,
        las = 2,
        ylab = "log2 intensity")

pheatmap(cor(data_matrix, method = "spearman"))

#plot PCA
mds <- plotMDS(
  data_matrix, 
  gene.selection = "common", 
  var.explained = TRUE
)

tibble(
  pc = 1:7,
  variance = mds$var.explained[1:7] * 100
) %>% 
  ggplot(aes(pc, variance)) +
  geom_col() +
  geom_text(aes(label = round(variance, digits = 2)), vjust = -0.2) +
  xlab("principal component number") +
  ylab("variance (%)")



# visualisation of DGE results --------------------------------------------

## heatmaps ---------------------------------------------------------------

plot_comparison <- function(logfc_col, p_col, file) {
  df_filtered <-
    df_data_unind_vs_treated %>%
    filter({{p_col}} < 0.05)
  
  ann_row <-
    df_filtered %>%
    column_to_rownames("accession") %>% 
    mutate(
      log2FC = if_else({{logfc_col}} > 0, "positive", "negative"),
      .keep = "none"
    )
  
  png(
    file,
    height = 100,
    width = 100,
    unit = "mm",
    res = 600
  )
  
  pheatmap(
    data_matrix[df_filtered$accession,],
    cluster_cols = TRUE,
    show_rownames = FALSE,
    annotation_colors = list(log2FC = c(negative = "blue", positive = "red")),
    annotation_row = ann_row,
    scale = "row"
  )
  
  dev.off()
}

dev.off()

# aci vs control
plot_comparison(
  coef_t3alwofvs_control,
  p_value_adj_t3alwofvs_control,
  "figures/heatmap_aci_vs_cntrl_600dpi.png"
)

# hp vs control
plot_comparison(
  coef_t2hpylvs_control,
  p_value_adj_t2hpylvs_control,
  "figures/heatmap_hp_vs_cntrl_600dpi.png"
)

# lps vs control
plot_comparison(
  coef_t1lp_svs_control,
  p_value_adj_t1lp_svs_control,
  "figures/heatmap_lps_vs_cntrl_600dpi.png"
)


## Figure 2B volcano plot -------------------------------------------------

ggplot(
  df_data_unind_vs_treated,
  aes(coef_t2hpylvs_control, -log10(p_value_adj_t2hpylvs_control))
) +
  geom_point() +
  xlim(-3, 4) +
  ylab("-log10(padj)") +
  xlab("log2 fold change")

