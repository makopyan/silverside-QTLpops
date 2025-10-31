library(tidyverse)
library(fdrtool)
library(data.table)

JIGA.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_JIGA_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
LRSC.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_LRSC_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
CTVA.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_CTVA_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
OINC.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_OINC_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
MCNC.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_MCNC_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
PANY.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_PANY_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
MCCT.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_MCCT_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
DXMA.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_DXMA_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
MDME.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_MDME_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
HXNS.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_HXNS_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
MBNS.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_MBNS_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))
MAQU.maf <- read_delim("/workdir/arne/results/snp_datasets/mme_MAQU_filtsnps_maf.mafs.gz",delim = "\t", col_names = T, col_types = list(chromo = col_character() ))

MAF.df <- cbind(JIGA.maf[,c(1,2,6)], 
                LRSC.maf[,6], 
                MCNC.maf[,6], 
                OINC.maf[,6], 
                CTVA.maf[,6],
                PANY.maf[,6],
                MCCT.maf[,6],
                DXMA.maf[,6],
                MDME.maf[,6],
                MBNS.maf[,6],
                HXNS.maf[,6],
                MAQU.maf[,6])

names(MAF.df) <- c("chromo","position","jiga.MAF"
                   ,"lrsc.MAF"
                   ,"mcnc.MAF"
                   ,"oinc.MAF"
                   ,"ctva.MAF"
                   ,"pany.MAF"
                   ,"mcct.MAF"
                   ,"dxma.MAF"
                   ,"mdme.MAF"
                   ,"mbns.MAF"
                   ,"hxns.MAF"
                   ,"maqu.MAF")
MAF.df <- as_tibble(MAF.df)

MAF.df.snps <- MAF.df[,c(1,2)]
MAF.df.freq <- MAF.df[,c(3:14)]

MAF.t <- t(MAF.df.freq)
MAF.t <- as_tibble(MAF.t)
MAF.t2 <- type_convert(MAF.t)

# Add column with latitude per population:
lat <- c(31.0, 33.9,34.7,35.8,37.9,40.8,41.3,42.0,44.4,45.3,44.6,47.4)
MAF.lat.df <- cbind(lat, MAF.t2)


# GLM function to loop over SNPs
proc_glm <- function(snps) {
  univariate <- glm(snps ~ MAF.lat.df$lat , family = gaussian, weights = NULL)
  return(coef(summary(univariate))[,4])
}

# run GLM outputting only p-value for lat
#glm_res <- as_tibble(t(as_tibble(lapply(MAF.lat.df[2:ncol(MAF.lat.df)], proc_glm))[2,]))
glm_res <- as_tibble(t(as_tibble(lapply(MAF.lat.df[2:10], proc_glm))[2,]))


# Perform BH correction
glm_res.padj <- as_tibble(p.adjust(glm_res$V1, method = "BH"))
names(glm_res.padj) <- "p.adj"

# Perform FDR correction
glm_res.fdr <- fdrtool(glm_res$V1, statistic="pvalue",
                        plot=F, color.figure=F, verbose=TRUE,
                        cutoff.method="fndr")


glm_res.fdr.pval <- as_tibble(glm_res.fdr$pval)
names(glm_res.fdr.pval) <- "pval"
glm_res.fdr.qval <- as_tibble(glm_res.fdr$qval)
names(glm_res.fdr.qval) <- "qval"
glm_res.fdr.lfdr <- as_tibble(glm_res.fdr$lfdr)
names(glm_res.fdr.lfdr) <- "lfdr"

glm.res.corr <- as_tibble(cbind(MAF.df.snps, glm_res.fdr.pval, glm_res.fdr.qval, glm_res.fdr.lfdr, glm_res.padj))

write_csv(glm.res.corr, "/workdir/arne/results/mme_assoc_results/glm_maf_lat/GLM_maf_vs_lat_bysnp_res.csv", col_names = T)