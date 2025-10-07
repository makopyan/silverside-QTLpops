# silverside-QTLpops
This repository contains code for analyzing the adaptive role of chromosomal inversions in Atlantic silverside populations. 

## QTL Analysis Files

Files for QTL analyses.

Location: [`QTL`](https://github.com/makopyan/silverside-QTLpops/tree/main/QTL)

## Wild Population Analysis Files

Files for population genomic analyses and inversion-trait associations.

Location: [`Wild`](https://github.com/makopyan/silverside-QTLpops/tree/main/Wild)

## Experimental Population Analysis Files

Files for generating and analyzing inversion karyotypes and their association with body size in artificial selection lines.

Location: [`Experiment`](https://github.com/makopyan/silverside-QTLpops/tree/main/Experiment)

- **[`all_bamlist_MixupsRemoved_expOnly`](https://github.com/makopyan/silverside-QTLpops/blob/main/Experiment/all_bamlist_MixupsRemoved_expOnly)** - List of sample BAM files
- **[`data_processing.md`](https://github.com/makopyan/silverside-QTLpops/blob/main/Experiment/data_processing.md)** - Documentation of data processing steps for genomic data
- **[`mme_exp_inv_local_pca.Rmd`](https://github.com/makopyan/silverside-QTLpops/blob/main/Experiment/mme_exp_inv_local_pca.Rmd)** - PCA script to infer inversion karyotypes 
- **[`mme_exp_main_linkage_inv_gts.csv`](https://github.com/makopyan/silverside-QTLpops/blob/main/Experiment/mme_exp_main_linkage_inv_gts.csv)** - Input file with inversion karyotypes and trait data
- **[`test_plot_exp_pxg.R`](https://github.com/makopyan/silverside-QTLpops/blob/main/Experiment/test_plot_exp_pxg.R)** - R script for testing and plotting phenotype-genotype associations
