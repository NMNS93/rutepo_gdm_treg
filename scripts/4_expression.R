# Run differential expression testing within clusters on GDM vs Controls

# Setup
source("code/differential_expression_v2.R")
source("code/stacked_violin.R")
source("code/plots.R")

library(qs)
library(fs)
library(dplyr)
library(readr)
library(stringr)
library(parallel)
library(data.table)

outdir <- here::here("gdm/results_v3")
figdir <- fs::dir_create("gdm/figures_v3/diff_expr")
gdm <- qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 <- gdm$cd4
treg <- gdm$treg
rm(gdm)
Idents(treg) <- "cluster"
Idents(cd4) <- "cluster"

# Differential expression testing ----

# Call differential expression
de_treg <- call_de_all(treg)
de_cd4 <- call_de_all(cd4)

# Add annotation for significant genes
de_treg <- left_join(de_treg, get_sig_de_genes_v3(de_treg))
de_cd4 <- left_join(de_cd4, get_sig_de_genes_v3(de_cd4))

# Save results
data.table::fwrite(de_treg, fs::path(outdir, "treg_within_cluster_de.csv.gz"))
data.table::fwrite(de_cd4, fs::path(outdir, "cd4_within_cluster_de.csv.gz"))
qs::qsave(list("wta" = de_treg, "wca" = de_cd4), fs::path(outdir, "we_treg_cln.qs"))
