# Run clustering and dimension reduction using best parameters identified in grid search

# Library ----
suppressPackageStartupMessages(source("code/preprocessing_qc.R"))
suppressPackageStartupMessages(source("code/plots.R"))
outdir <- fs::dir_create(here::here("gdm/results"))
cache_dir <- fs::dir_create(fs::path(outdir, "cache"))
figdir <- fs::dir_create(here::here("gdm/figures/clustering"))

# Variables ----
gdm_seurat <- readRDS(fs::path(outdir, "gdm_seurat.rds"))
treg <- gdm_seurat$treg
cd4 <- gdm_seurat$cd4
rm(gdm_seurat)

# Main ----

# Filter clusters with fewer than 20 cells per patient (median) and cluster
treg <- filter_cells_pat_median(treg, 20)
cd4 <- filter_cells_pat_median(cd4, 20)
treg <- RunDimRed(RunSeuratClustering(treg, resolution = 0.4), components = 3, neighbors = 20, perplexity = 15)
cd4 <- RunDimRed(RunSeuratClustering(cd4, resolution = 0.9), components = 2, neighbors = 7, perplexity = 35)

# Show number of cells per patient and median per cluster
p_cpp_med <- plot_med_per_pat_clust(cells_per_pat_clust(treg, "treg")) /
  plot_med_per_pat_clust(cells_per_pat_clust(cd4, "cd4"))
ggsave(fs::path(figdir, "cells_per_pat_median.png"), p_cpp_med, width = 20, height = 10)
ggsave(fs::path(figdir, "cells_per_pat_raw.png"),
  show_count_pat_clust(treg) / show_count_pat_clust(cd4),
  width = 12, height = 20
)

# Show DimPlots across all four clusters after filtering and reclustering
all_dimplots_filt <- show_dimplots(treg, cd4)
ggsave(fs::path(figdir, "all_dimplots_filt.png"), all_dimplots, width = 10, height = 10)

# Save data objects
qs::qsave(list("treg" = treg, "cd4" = cd4), fs::path(outdir, "gdm_seurat_reclust_filter.qs"))
