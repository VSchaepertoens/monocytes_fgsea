## ---------------------------
##
## Script name: analyse_dge
##
## Purpose of script: Using limma package, analyse differentially expressed proteins in 
## treated (LPS, Hp, Aci) and uninduced monocytes
##
## Author: Dr. Christof Regl & Dr. Veronika Schäpertöns
##
## Date Created: 16.03.2023
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
#LIMMA-------------------------------------------------------------- 


library(limma)

#import data (already log2 transformed normalized and batch corrected by limma in Perseus)-------


df <- read.delim("data/proteinGroups_log2_substractMedian_Batchcorrected.csv", 
                 sep = ';', 
                 header = TRUE)
dim(df)
df

#design matrix------------------------------------------------------------

sample <- as.factor(c("Control", "Control", "Control", "T1LPS", "T1LPS", "T1LPS", "T2Hpyl", "T2Hpyl", "T2Hpyl", "T3Alwof", "T3Alwof", "T3Alwof"))
replicate <- as.factor(rep(c("1","2","3", "1", "2", "3", "1", "2", "3", "1","2","3") ))


design.matrix <- model.matrix(~0+sample+replicate)
colnames(design.matrix) <- gsub("sample", 
                                "", 
                                colnames(design.matrix))
design.matrix




# Contrasts for LIMMA-----------------------------------------------------

contr.matrix <- makeContrasts(
  T1LPSvsControl = T1LPS  - Control,
  T2HpylvsControl = T2Hpyl - Control,
  T3AlwofvsControl = T3Alwof - Control,
  T2HpylvsT1LPS = T2Hpyl - T1LPS, 
  T3AlwofvsT1LPS = T3Alwof - T1LPS, 
  T3AlwofvsT2Hpyl = T3Alwof - T2Hpyl, 
  levels = colnames(design.matrix))
contr.matrix


# LIMMA---------------------------------------------------------------------
fit1 <- lmFit(df, design.matrix)
fit2 <- contrasts.fit(fit1, contrasts = contr.matrix)
fit3 <- eBayes(fit2)

summary(decideTests(fit3))

topTable(fit3)

results <- data.frame(fit3)

write.fit(fit3, 
          adjust = "BH",
          sep = ";", 
          row.names = FALSE, 
          file = "analysis/proteinGroups_log2_substractMedian_Batchcorrected_Limma.csv")
