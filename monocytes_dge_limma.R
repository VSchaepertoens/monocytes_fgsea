# differential expression analysis using limma

library(tidyverse)
library(janitor)
library(pheatmap)
library(limma)


# load dataset ------------------------------------------------------------

df_log2_all <- read_csv("data/log2_substractMed_all_replicate.csv",
                        col_names = TRUE) %>%
  column_to_rownames(var = "...1") %>%
  clean_names()

data_matrix <- as.matrix(df_log2_all)  


# create meta data table --------------------------------------------------

design_matrix <- data.frame(sample_name = colnames(df_log2_all)) %>% 
  separate(
    sample_name, 
    into = c("treatment", "replicate"), 
    sep = "_", 
    remove = FALSE)


# data exploration PCA & correlations ---------------------------------------

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

# create design matrix ----------------------------------------------------

des <- copy(design_matrix)
unique(des$treatment)

des$treatment <- factor(des$treatment,
                        levels = c("uninduced", "lps", "hpyl", "alwof"))

des <- model.matrix(~treatment,
                    data = des)

# model fit ---------------------------------------------------------------

fit <- lmFit(data_matrix, des)
fit <- eBayes(fit)

# get results -------------------------------------------------------------

coefs <- grep("treatment", colnames(coef(fit)), value = TRUE)
res <- data.table()
for (coefx in coefs) {
  res <- rbind(res, data.table(
    topTable(fit, coef = coefx, number = nrow(data_matrix)), 
    keep.rownames = TRUE,
    coef = gsub("treatment", "", coefx)
  ))
}
res[coef == "hpyl"][adj.P.Val < 0.05]
res[coef == "alwof"][adj.P.Val < 0.05]
res[coef == "lps"][adj.P.Val < 0.05]









