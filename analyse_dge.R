##
## Script name: analyse_dge
##
## Purpose of script: Analyse differentially expressed proteins in  treated (LPS, Hp, Aci)
## and uninduced monocytes via the limma package
##
## Author: Dr. Christof Regl & Dr. Veronika Schäpertöns
##
## Date Created: 16.03.2023
##
## Copyright (c) Veronika Schäpertöns, 2023
## Email: veronika.schaepertoens@plus.ac.at
##
##
## Notes:
##   
##

library(limma)
library(fs)



# Import data -------------------------------------------------------------

# data was already log2 transformed, normalized and batch corrected
# by limma in Perseus

df <- read.delim("data/proteinGroups_log2_substractMedian_Batchcorrected.csv", 
                 sep = ';', 
                 header = TRUE)
dim(df)
df



# Create design matrix ----------------------------------------------------

sample <- as.factor(rep(c("Control", "T1LPS", "T2Hpyl", "T3Alwof"), each = 3))
replicate <- as.factor(rep(1:3, 4))


design.matrix <- model.matrix(~0 + sample + replicate)
colnames(design.matrix) <- gsub("sample", 
                                "", 
                                colnames(design.matrix))
design.matrix



# Contrasts for LIMMA -----------------------------------------------------

contr.matrix <- makeContrasts(
  T1LPSvsControl = T1LPS  - Control,
  T2HpylvsControl = T2Hpyl - Control,
  T3AlwofvsControl = T3Alwof - Control,
  T2HpylvsT1LPS = T2Hpyl - T1LPS, 
  T3AlwofvsT1LPS = T3Alwof - T1LPS, 
  T3AlwofvsT2Hpyl = T3Alwof - T2Hpyl, 
  levels = colnames(design.matrix)
)
contr.matrix



# Run LIMMA ---------------------------------------------------------------

fit1 <- lmFit(df, design.matrix)
fit2 <- contrasts.fit(fit1, contrasts = contr.matrix)
fit3 <- eBayes(fit2)

summary(decideTests(fit3))

topTable(fit3)

results <- data.frame(fit3)



# Save results ------------------------------------------------------------

dir_create("analysis")

write.fit(fit3, 
          adjust = "BH",
          sep = ",", 
          row.names = FALSE, 
          file = "analysis/results_limma.csv")
