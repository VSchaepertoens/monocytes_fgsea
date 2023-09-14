# differential expression analysis using limma

library(tidyverse)
library(janitor)
library(pheatmap)
library(limma)
library(data.table)
library(readxl)


# exploration of data -----------------------------------------------------


## load dataset ------------------------------------------------------------

df_log2_all <- read_xlsx("analysis/Monocytes_Experiment_Results_20230316.xlsx", 
                           sheet = "log2_substractMed_all_replicate") %>%
  column_to_rownames(var = "Accession") %>%
  clean_names()

data_matrix <- as.matrix(df_log2_all)  

## create meta data table --------------------------------------------------

design_matrix <- data.frame(sample_name = colnames(df_log2_all),
                            label = rep(c("Uninduced","LPS","Hp","Aci"), each = 3)
                            ) %>% 
  separate(
    sample_name, 
    into = c("treatment", "replicate"), 
    sep = "_", 
    remove = FALSE)


## data exploration PCA & correlations ---------------------------------------

plot(density(data_matrix))

boxplot(data_matrix,
        las = 2,
        names = design_matrix$sample_name,
        ylab = "log2 intensity")

pheatmap(cor(data_matrix,method = "spearman"))

#plot PCA
mds <- plotMDS(data_matrix, 
        labels = design_matrix$sample_name, 
        gene.selection = "common", 
        var.explained = TRUE)

var_explained <- as.data.frame(mds$var.explained[1:7]*100)
colnames(var_explained) <- c("variance")
ggplot(var_explained, aes(x = rownames(var_explained), y = variance)) +
  geom_col() +
  geom_text(aes(label = round(variance, digits = 2)),vjust = -0.2) +
  xlab("principal component number") +
  ylab("variance (%)")

# # create design matrix ----------------------------------------------------
# 
# des <- copy(design_matrix)
# unique(des$treatment)
# 
# des$treatment <- factor(des$treatment,
#                         levels = c("uninduced", "lps", "hpyl", "alwof"))
# 
# des <- model.matrix(~treatment,
#                     data = des)
# 
# # model fit ---------------------------------------------------------------
# 
# fit <- lmFit(data_matrix, des)
# fit <- eBayes(fit)
# 
# # get results -------------------------------------------------------------
# 
# coefs <- grep("treatment", colnames(coef(fit)), value = TRUE)
# res <- data.table()
# for (coefx in coefs) {
#   res <- rbind(res, data.table(
#     topTable(fit, coef = coefx, number = nrow(data_matrix)), 
#     keep.rownames = TRUE,
#     coef = gsub("treatment", "", coefx)
#   ))
# }
# res[coef == "hpyl"][adj.P.Val < 0.05]
# res[coef == "alwof"][adj.P.Val < 0.05]
# res[coef == "lps"][adj.P.Val < 0.05]





# visualisation of DGE results --------------------------------------------


## load tables with sig proteins from analyse_dge --------------------------

df_data_unind_vs_treated <- read_xlsx("analysis/Monocytes_Experiment_Results_20230316.xlsx", 
                                      sheet = "log2_fold_changes_all_treatment") %>%
  clean_names()


## lps_vs_cntrl ------------------------------------------------------------

filtered <- df_data_unind_vs_treated %>%
  filter(p_value_adj_t1lp_svs_control < 0.05) %>%
  as.data.table()

ann.row <- with(filtered, data.frame(row.names = accession, log2FC = ifelse(coef_t1lp_svs_control > 0, 1, -1)))

ann_col = data.frame(label_sample = paste(design_matrix$label, design_matrix$replicate, sep = "_")) 

png("figures/heatmap_lps_vs_cntrl_600dpi.png",
    height = 100,
    width = 100,
    unit = "mm",
    res = 600)

pheatmap(data_matrix[filtered$accession,],
         cluster_cols = TRUE,
         show_rownames = FALSE,
         labels_col = ann_col$label_sample,
         annotation_colors = list(log2FC = c("blue","red")),
         annotation_row = ann.row,
         scale = "row",
         width = ,
         height = )

dev.off()

## hp_vs_cntrl -------------------------------------------------------------
### Figure 2C heatmap -------------------------------------------------------

filtered <- df_data_unind_vs_treated %>%
  filter(p_value_adj_t2hpylvs_control < 0.05) %>%
  as.data.table()

ann.row <- with(filtered, data.frame(row.names = accession, log2FC = ifelse(coef_t2hpylvs_control > 0, 1, -1)))

ann_col = data.frame(label_sample = paste(design_matrix$label, design_matrix$replicate, sep = "_")) 

png("figures/heatmap_hp_vs_cntrl_600dpi.png",
    height = 100,
    width = 100,
    unit = "mm",
    res = 600)

pheatmap(data_matrix[filtered$accession,],
         cluster_cols = TRUE,
         show_rownames = FALSE,
         labels_col = ann_col$label_sample,
         annotation_colors = list(log2FC = c("blue","red")),
         annotation_row = ann.row,
         scale = "row",
         width = ,
         height = )
dev.off()


### Figure 2B volcano plot --------------------------------------------------

ggplot(df_data_unind_vs_treated, aes(x = coef_t2hpylvs_control, y = -log10(p_value_adj_t2hpylvs_control))) +
  geom_point() +
  xlim(-3, 4) +
  ylab("-log10(padj)") +
  xlab("log2 fold change")


## aci_vs_cntrl ------------------------------------------------------------

filtered <- df_data_unind_vs_treated %>%
  filter(p_value_adj_t3alwofvs_control < 0.05) %>%
  as.data.table()

ann.row <- with(filtered, data.frame(row.names = accession, log2FC = ifelse(coef_t3alwofvs_control > 0, 1, -1)))

ann_col = data.frame(label_sample = paste(design_matrix$label, design_matrix$replicate, sep = "_")) 

png("figures/heatmap_aci_vs_cntrl_600dpi.png",
    height = 100,
    width = 100,
    unit = "mm",
    res = 600)

pheatmap(data_matrix[filtered$accession,],
         cluster_cols = TRUE,
         show_rownames = FALSE,
         labels_col = ann_col$label_sample,
         annotation_colors = list(log2FC = c("blue","red")),
         annotation_row = ann.row,
         scale = "row",
         width = ,
         height = )
dev.off()
