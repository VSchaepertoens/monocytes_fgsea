# monocytes_fgsea
functional class scoring of monocyte data

## Download data

Create a folder `data` that will contain raw data downloaded from the following link: [link to data]

## Main workflow

Run these R scripts in the given in order to generate all files & plots

-   [download_enrichr_databases.R](download_enrichr_databases.R) - Downloads respective databases from Enrichr for fGSEA analysis
-   [monocytes_gsea.R](monocytes_gsea.R) - Runs fGSEA analysis on all data and plots figures
