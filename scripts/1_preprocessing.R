# Preprocess raw data for RUTEPO scRNA analysis.

# Library ----
suppressPackageStartupMessages(source("code/preprocessing_qc.R"))
suppressPackageStartupMessages(source("code/plots.R"))
suppressPackageStartupMessages(source("code/differential_expression.R"))
suppressPackageStartupMessages(source("code/cluster_umap_tsne_search.R"))

# Variables ----
outdir <- fs::dir_create(here::here("gdm/results"))
cache_dir <- fs::dir_create(fs::path(outdir, "cache"))
figdir <- fs::dir_create(here::here("gdm/figures"))
cellranger_csv <- here::here("data/cellranger/220802_Shangaris_CellRangerDirs.csv")
hk_csv <- here::here("data/seurat/ref/Housekeeping_GenesHuman.csv")
meta_csv <- here::here("data/metadata/220721_metadata_cln.csv")

# Main ----

log_info("Starting preprocessing")

# Run preprocessing. Returns an object with $treg and $cd4 slots
gdm_seurat <- RunGDMPreproccessingQC(cellranger_csv, hk_csv, meta_csv, cache_dir)
saveRDS(gdm_seurat, fs::path(outdir, "gdm_seurat.rds"))

log_info("preprocessing complete")
