# Perform sensitivity analysis excluding L215-A-AHH05

# Preprocessing ----

# Setup
library(dplyr)
library(Seurat)
library(fs)
library(qs)

# Input
f_seurat <- "gdm/results/gdm_seurat.rds"

# Output
f_seurat_out <- "gdm/results/sensitivity/gdm_seurat.qs"

# Create output directory if it doesn't exist
fs::dir_create(
  fs::path_dir(f_seurat_out)
)

# Custom function to filter Seurat object
filter_ <- function(x) {
  cells <- WhichCells(x, expression = multi_q != "L215-A-AHH05")
  subset(x, cells = cells)
}

# Main
seu <- readRDS(f_seurat)
seu$treg <- filter_(seu$treg)
seu$cd4 <- filter_(seu$cd4)
qs::qsave(seu, f_seurat_out)

# Clustering and dimensionality reduction ----

# Repeat analysis using parameters identified in 1a

# Setup
library(qs)
library(Matrix)
library(Seurat)
suppressPackageStartupMessages(source("code/preprocessing_qc.R"))
suppressPackageStartupMessages(source("code/plots.R"))

outdir <- fs::dir_create(here::here("gdm/results/sensitivity"))
cache_dir <- fs::dir_create(fs::path(outdir, "cache"))
figdir <- fs::dir_create(here::here("gdm/figures/sensitivity/clustering"))

# Load datasets
gdm_seurat <- qs::qread(fs::path(outdir, "gdm_seurat.qs"))
treg <- gdm_seurat$treg
cd4 <- gdm_seurat$cd4
rm(gdm_seurat)

# Recluster with chosen parameters
treg <- RunDimRed(RunSeuratClustering(treg, resolution = 0.4), components = 3, neighbors = 20, perplexity = 15)
cd4 <- RunDimRed(RunSeuratClustering(cd4, resolution = 0.9), components = 2, neighbors = 7, perplexity = 35)

# Show DimPlots across all four clusters
all_dimplots <- show_dimplots(treg, cd4)
ggsave(fs::path(figdir, "all_dimplots.png"), all_dimplots, width = 10, height = 10)

# Drop clusters with fewer than 20 cells per patient (median) and repeat DimReg
treg <- filter_cells_pat_median(treg, 20)
cd4 <- filter_cells_pat_median(cd4, 20)
treg <- RunDimRed(RunSeuratClustering(treg, resolution = 0.4), components = 3, neighbors = 20, perplexity = 15)
cd4 <- RunDimRed(RunSeuratClustering(cd4, resolution = 0.9), components = 2, neighbors = 7, perplexity = 35)

# Show number of cells per patient and median per cluster
p_cpp_med <- plot_cells_per_pat_median(treg, cd4) # Median plots
ggsave(fs::path(figdir, "cells_per_pat_median.png"), p_cpp_med, width = 20, height = 10)
ggsave(fs::path(figdir, "cells_per_pat_raw.png"), # Bar plots
  show_count_pat_clust(treg) / show_count_pat_clust(cd4),
  width = 12, height = 20
)

# Show DimPlots across all four clusters after filtering and reclustering
all_dimplots_filt <- show_dimplots(treg, cd4)
ggsave(fs::path(figdir, "all_dimplots_filt.png"), all_dimplots, width = 10, height = 10)

# Save new objects
qs::qsave(list("treg" = treg, "cd4" = cd4), fs::path(outdir, "gdm_seurat_reclust_filter.qs"))
log_info("Analysis complete.")

# Subtype tracking ----
# Calculate the percentage of cells assigned to the same subtype after dropping L215-A-AHH05

source("code/helpers.R")
source("code/plots.R")
source("code/subtyping_de.R")
library(dplyr)
library(data.table)
library(ggplot2)

# Define a function that counts the cell re-assignments for each group
sens_count <- function(orig, obj, str_cluster) {
  # Get name
  cluster_name <- orig@meta.data %>%
    dplyr::filter(seurat_clusters == str_cluster) %>%
    dplyr::slice(1) %>%
    dplyr::pull(cluster)

  # Find cells
  cells <- WhichCells(orig, expression = seurat_clusters == str_cluster)

  # Select cells in new dataset
  sub <- subset(obj, cells = cells)

  # Count assignments
  n_assign <- dplyr::count(sub@meta.data, seurat_clusters) %>%
    dplyr::rename(new_cluster = seurat_clusters) %>%
    dplyr::mutate(
      old_cluster = str_cluster,
      old_cluster_name = cluster_name,
      percentages = 100 * n / sum(n)
    )

  return(n_assign)
}

# Args
origdir <- fs::path("gdm/results/")
outdir <- fs::dir_create("gdm/results/sensitivity/")
figdir <- fs::dir_create("gdm/figures/sensitivity/subtyping")

# Load dataset
gdm <- qs::qread(fs::path(outdir, "gdm_seurat_reclust_filter.qs"))
treg <- gdm$treg
cd4 <- gdm$cd4
rm(gdm)
DefaultAssay(treg) <- "SCT"
DefaultAssay(cd4) <- "SCT"

# Load original dataset
ogdm <- qs::qread(fs::path(origdir, "gdm_subtyped.qs"))
otreg <- ogdm$treg
ocd4 <- ogdm$cd4
rm(ogdm)
DefaultAssay(otreg) <- "SCT"
DefaultAssay(ocd4) <- "SCT"

# Apply counts
treg_clusters <- (unique(otreg$seurat_clusters))
treg_re <- data.table::rbindlist(
  lapply(
    treg_clusters,
    sens_count,
    orig = otreg,
    obj = treg
  )
)
treg_re$source <- "treg"

cd4_clusters <- (unique(ocd4$seurat_clusters))
cd4_re <- data.table::rbindlist(
  lapply(
    cd4_clusters,
    sens_count,
    orig = ocd4,
    obj = cd4
  )
)
cd4_re$source <- "cd4"

# Save datasets
res <- rbind(treg_re, cd4_re)
data.table::fwrite(res, fs::path(outdir, "subtype_sensitivity_scores.csv"))

# Visualise the scores in the new clusters
rmax <- res %>%
  dplyr::group_by(old_cluster, old_cluster_name, source) %>%
  dplyr::slice_max(order_by = percentages, n = 1)
rmax

library(patchwork)

p1 <- ggplot(dplyr::filter(rmax, source == "treg")) +
  geom_col(aes(x = old_cluster_name, y = percentages, fill = new_cluster)) +
  theme_bw() +
  scale_fill_viridis_d() +
  labs(x = "", y = "Cells retained in cluster (%)", fill = "Sensitivity\nAnalysis\nCluster", title = "Treg")

p2 <- ggplot(dplyr::filter(rmax, source == "cd4")) +
  geom_col(aes(x = old_cluster_name, y = percentages, fill = new_cluster)) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "", y = "Cells retained in cluster (%)", fill = "Sensitivity\nAnalysis\nCluster", title = "CD4")
ggsave(fs::path(figdir, "sensitivity_analysis_group_percentages.png"), p1 / p2, width = 8, height = 8)

# Subtyping ----

# Perform differential expression testing between clusters to determine which genes identify each cluster.
source("code/helpers.R")
source("code/plots.R")
source("code/subtyping_de.R")
library(dplyr)

# Args
origdir <- fs::path("gdm/results/")
outdir <- fs::dir_create("gdm/results/sensitivity/")
figdir <- fs::dir_create("gdm/figures/sensitivity/subtyping")

# Load dataset
gdm <- qs::qread(fs::path(outdir, "gdm_seurat_reclust_filter.qs"))
treg <- gdm$treg
cd4 <- gdm$cd4
rm(gdm)
DefaultAssay(treg) <- "SCT"
DefaultAssay(cd4) <- "SCT"

# Load original dataset
ogdm <- qs::qread(fs::path(origdir, "gdm_subtyped.qs"))
otreg <- ogdm$treg
ocd4 <- ogdm$cd4
rm(ogdm)
DefaultAssay(otreg) <- "SCT"
DefaultAssay(ocd4) <- "SCT"

# Get meta
tmeta <-
  # Set clusters
  treg_cluster <- c("Naive-1", "Naive-2", "Effector-1", "Effector-2", "Effector-3", "Humanin+", "MALAT1+")
names(treg_cluster) <- seq(0, 6)
treg_clust_group <- c("Naive", "Naive", "Effector", "Effector", "Effector", "Undefined", "Undefined")
names(treg_clust_group) <- seq(0, 6)
cd4_cluster <- c("Naive-2", "Naive-1", "Memory-1", "Naive-3", "Memory-2", "Naive-4", "Effector-1", "Humanin+", "FOXP3+")
names(cd4_cluster) <- seq(0, 8)
cd4_clust_group <- c("Naive", "Naive", "Memory", "Naive", "Memory", "Naive", "Effector", "Humanin+", "FOXP3+")
names(cd4_clust_group) <- seq(0, 8)

# Add metadata
treg <- AddMetaData(treg, treg_cluster[Idents(treg)], col.name = "cluster")
treg <- AddMetaData(treg, treg_clust_group[Idents(treg)], col.name = "clust_group")
cd4 <- AddMetaData(cd4, cd4_cluster[Idents(cd4)], col.name = "cluster")
cd4 <- AddMetaData(cd4, cd4_clust_group[Idents(cd4)], col.name = "clust_group")

# Save
gdm_subtyped <- list(treg = treg, cd4 = cd4)
qs::qsave(gdm_subtyped, fs::path(outdir, "gdm_subtyped.qs"))

# Differential expression ----

# Run differential expression testing within clusters on GDM vs Controls

# Setup
source("code/differential_expression.R")
source("code/plots.R")

library(qs)
library(fs)
library(dplyr)
library(readr)
library(stringr)
library(parallel)
library(data.table)

outdir <- here::here("gdm/results/sensitivity")
figdir <- fs::dir_create("gdm/figures/sensitivity/diff_expr")
gdm <- qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 <- gdm$cd4
treg <- gdm$treg
rm(gdm)
Idents(treg) <- "cluster"
Idents(cd4) <- "cluster"

# Differential expression testing ----
# Set label field
treg$label <- "treg"
cd4$label <- "cd4"

# Call differential expression
de_treg <- call_de_all(treg)
de_cd4 <- call_de_all(cd4)
qs::qsave(de_treg, fs::path(outdir, "de_treg_list.qs"))
qs::qsave(de_cd4, fs::path(outdir, "de_cd4_list.qs"))

# Add annotation for significant genes
de_treg <- left_join(de_treg, get_sig_de_genes_v3(de_treg))
de_cd4 <- left_join(de_cd4, get_sig_de_genes_v3(de_cd4))

# Save results
data.table::fwrite(de_treg, fs::path(outdir, "treg_within_cluster_de.csv.gz"))
data.table::fwrite(de_cd4, fs::path(outdir, "cd4_within_cluster_de.csv.gz"))
qs::qsave(list("wta" = de_treg, "wca" = de_cd4), fs::path(outdir, "we_treg_cln.qs"))

# Expression matching ----

# For each cluster, count the number of significant genes that remain significantly DE after sensitivity analysis.
library(qs)
library(dplyr)

outdir <- here::here("gdm/results/sensitivity")
prevdir <- here::here("gdm/results/")
figdir <- here::here("gdm/figures/sensitivity")

de <- qs::qread(fs::path(outdir, "we_treg_cln.qs"))
ode <- qs::qread(fs::path(prevdir, "we_treg_cln.qs"))

# For each gene used in the final model, determine whether it is present in the DE genes.
mtreg <- merge(
  de$wta %>% dplyr::select(gene, label, cluster, is_signif),
  ode$wta %>% dplyr::select(gene, label, cluster, was_signif = is_signif),
  all = T
) %>% dplyr::filter(was_signif == T)
streg <- mtreg %>%
  dplyr::group_by(label, cluster) %>%
  dplyr::summarise(ssignif = sum(is_signif, na.rm = T), swasignif = sum(was_signif, na.rm = T)) %>%
  tidyr::pivot_longer(cols = c("ssignif", "swasignif"))
streg

mcd4 <- merge(
  de$wca %>% dplyr::select(gene, label, cluster, is_signif),
  ode$wca %>% dplyr::select(gene, label, cluster, was_signif = is_signif),
  all = T
) %>% dplyr::filter(was_signif == T)
scd4 <- mcd4 %>%
  dplyr::group_by(label, cluster) %>%
  dplyr::summarise(ssignif = sum(is_signif, na.rm = T), swasignif = sum(was_signif, na.rm = T)) %>%
  tidyr::pivot_longer(cols = c("ssignif", "swasignif"))
scd4

library(ggplot2)
library(patchwork)
p1 <- ggplot(streg, aes(x = cluster, y = value, fill = name)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("grey", "blue"), labels = c("ssignif" = "Sensitivity\nAnalysis", "swasignif" = "Original", title = "")) +
  labs(x = "", y = "Significant GDM DEGs\nin original analysis", title = "Treg", fill = "")
p2 <- ggplot(scd4 %>% dplyr::filter(cluster != "Effector-2"), aes(x = cluster, y = value, fill = name)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 9), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("grey", "red"), labels = c("ssignif" = "Sensitivity\nAnalysis", "swasignif" = "Original", title = "")) +
  labs(x = "", y = "Significant GDM DEGs\nin original analysis", title = "CD4", fill = "") # False CD4 E2 result in ODE
p1 / p2
ggsave(fs::path(figdir, "sensitivity_de.png"), p1 / p2, width = 8.5, height = 6)

# Finally, plot the log2Fold change for DEGs that appear in both
ltreg <- merge(
  de$wta %>% dplyr::filter(is_signif == T) %>% dplyr::select(avg_log2FC, gene, label, cluster),
  ode$wta %>% dplyr::filter(is_signif == T) %>% dplyr::select(orig_log2FC = avg_log2FC, gene, label, cluster),
)
lcd4 <- merge(
  de$wca %>% dplyr::filter(is_signif == T) %>% dplyr::select(avg_log2FC, gene, label, cluster),
  ode$wca %>% dplyr::filter(is_signif == T) %>% dplyr::select(orig_log2FC = avg_log2FC, gene, label, cluster),
)

l1 <- ggplot(ltreg, aes(x = orig_log2FC, y = avg_log2FC)) +
  geom_point(alpha = 0.5, colour = "blue") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Log2FC (Original)", y = "Log2FC (Sensitivity Analysis)", title = "Treg")
l2 <- ggplot(lcd4, aes(x = orig_log2FC, y = avg_log2FC)) +
  geom_point(alpha = 0.5, colour = "red") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Log2FC (Original)", y = "Log2FC (Sensitivity Analysis)", title = "CD4")
ggsave(fs::path(figdir, "sensitivity_log2fc.png"), l1 + l2, width = 8, height = 4)
