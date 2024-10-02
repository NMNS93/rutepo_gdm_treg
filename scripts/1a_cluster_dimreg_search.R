# Explore parameters for clustering and dimension reduction

suppressPackageStartupMessages(source("code/plots.R"))
suppressPackageStartupMessages(source("code/differential_expression.R"))
suppressPackageStartupMessages(source("code/cluster_umap_tsne_search.R"))

outdir <- fs::dir_create(here::here("gdm/results_v3"))
figdir <- fs::dir_create(here::here("gdm/figures_v3"))
gdm_seurat <- readRDS(fs::path(outdir, "gdm_seurat.rds"))
treg <- gdm_seurat$treg
treg <- FindClusters(treg, resolution = 0.4) # Default high clust
cd4 <- gdm_seurat$cd4
cd4 <- FindClusters(cd4, resolution = 0.4) # Default high clust

# Optimise clustering and dimension reduction parameters
log_info("Starting grid search")
reclust_dir <- fs::dir_create(fs::path(figdir, "grid_search"))
search_umap(treg, "treg", reclust_dir, cores = 12)
search_clusters(treg, "treg", reclust_dir, cores = 12)
search_tsne(treg, "treg", reclust_dir, cores = 12)
bind_figs(reclust_dir, "treg_umap")
# Grid search CD4 params
search_umap(cd4, "cd4", reclust_dir, cores = 12)
search_clusters(cd4, "cd4", reclust_dir, cores = 12)
search_tsne(cd4, "cd4", reclust_dir, cores = 12)
bind_figs(reclust_dir, "cd4")
log_info("Grid search complete")
