# Monocytes_fgsea
Differential expression analysis and functional class scoring of monocyte data

## Download data

Create a folder `data` that will contain raw data downloaded from the following link: [link to data]

## Main workflow

Run these R scripts in the given in order to generate all files & plots

-   [analyse_dge.R](analyse_dge.R) - Loads log2 transformed, normalized, & batch-corrected dataset and analyses differential expression of proteins in treated (LPS, Hp, Aci) and uninduced monocytes
-   [monocytes_exploration_dge_heatmaps.R](monocytes_exploration_dge_heatmaps.R) - Loads log2 transformed, normalized, & batch-corrected dataset and perfoms explorative data analysis (PCA, Pearson correlations), loads results of DGE analysis and plots heatmaps of raw expressions only for significantly regulated proteins in lps_vs_uninduced, hp_vs_uninduced, aci_vs_uninduced conditions
-   [download_enrichr_databases.R](download_enrichr_databases.R) - Downloads respective databases from Enrichr for fGSEA analysis
-   [monocytes_gsea.R](monocytes_gsea.R) - Runs fGSEA analysis on all data and plots figures
