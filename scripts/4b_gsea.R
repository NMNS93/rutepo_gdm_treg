# Run GSEA on hallmarks for all scRNA clusters with DE
#   Use both p_values and log2FC as ranking

source("code/gsea.R")
source("code/plots.R")
library(qs)
library(fs)
library(dplyr)
library(readr)
library(logger)
library(parallel)
library(gridExtra)

# Load data
outdir = here::here("gdm/results_v3")
figdir = fs::dir_create("gdm/figures_v3/gsea")
gdm = qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 = gdm$cd4; treg = gdm$treg; rm(gdm)

# Load differential expression
we <- qs::qread(fs::path(outdir,"we_treg_cln.qs"))
we_treg = we$wta; we_cd4 = we$wca

# Load hallmark GMT files for initial GSEA
gmts <- fs::dir_ls("data/gsea_gene_sets/", glob="*h.all.v7.5.1.symbols.gmt")

# Run GSEA for each cluster
log_info("Treg gsea")
we_treg_pass <- dplyr::filter(we_treg, pct.1>=0.05 & pct.2 >=0.05)
tga <- gsea_runner_v3(we_treg_pass, gmts)
qs::qsave(tga, fs::path(outdir, "gsea_treg.qs"))

log_info("CD4 gsea")
we_cd4_pass <- dplyr::filter(we_cd4, pct.1>=0.05 & pct.2 >=0.05)
cga <- gsea_runner_v3(we_cd4, gmts)
qs::qsave(cga, fs::path(outdir, "gsea_cd4.qs"))

