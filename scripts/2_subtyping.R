# Perform differential expression testing between clusters to determine which genes identify each cluster.
source("code/helpers.R")
source("code/plots.R")
source("code/subtyping_de.R")
library(dplyr)

# Setup
outdir <- fs::dir_create("gdm/results_v3/")
figdir <- fs::dir_create("gdm/figures_v3/subtyping")
gdm <- qs::qread(fs::path(outdir, "gdm_seurat_reclust_filter.qs"))
treg <- gdm$treg
cd4 <- gdm$cd4
rm(gdm)
DefaultAssay(treg) <- "SCT"
DefaultAssay(cd4) <- "SCT"

# Prepare samples for DE testing
treg <- PrepSCTFindMarkers(treg)
cd4 <- PrepSCTFindMarkers(cd4)

# Call subtyping DE
treg_markers <- subtyping_de(treg, "treg")
cd4_markers <- subtyping_de(cd4, "cd4")
subtyping_de_markers <- rbind(treg_markers, cd4_markers)
readr::write_csv(subtyping_de_markers, fs::path(outdir, "subtyping_de_markers.csv"))

# Filter and save top-performing markers
top_markers <- filter_for_top_markers(subtyping_de_markers)
readr::write_csv(top_markers, fs::path(outdir, "top_subtyping_de_markers.csv"))
qs::qsave(list("treg" = treg, "cd4" = cd4), fs::path(outdir, "gdm_seurat_sct_markers.qs"))

# Show the top 3 DE plots for each cluster (treg)
gdm <- qs::qread(fs::path(outdir, "gdm_seurat_sct_markers.qs"))
treg <- gdm$treg
cd4 <- gdm$cd4
rm(gdm)
de_treg <- readr::read_csv(fs::path(outdir, "top_subtyping_de_markers.csv")) %>% dplyr::filter(label == "treg")
p_de_treg_top3 <- ShowSCTClusterDeTop(treg, de_treg)
treg_de_idents <- seq_along(unique(de_treg$ident))
for (i in treg_de_idents) {
    ggsave(fs::path(figdir, paste0("p_treg_subtype_top3_", i, ".png")), p_de_treg_top3[[i]], width = 15, height = 6)
}

# Show the top 3 DE plots for each cluster (cd4)
de_cd4 <- readr::read_csv(fs::path(outdir, "top_subtyping_de_markers.csv")) %>% dplyr::filter(label == "cd4")
p_de_cd4_top3 <- ShowSCTClusterDeTop(cd4, de_cd4)
cd4_de_idents <- seq_along(unique(de_cd4$ident))
for (i in cd4_de_idents) {
    ggsave(fs::path(figdir, paste0("p_cd4_subtype_top3_", i, ".png")), p_de_cd4_top3[[i]], width = 15, height = 6)
}

# Markers
genes <- c(
    "CCR7", "CD28", "CD95", "CD45RA", "CCR4", "CCR6", "CXCR3", "CXCR5", "L-selectin",
    "CTLA-4", "ICOS", "Ki67", "BCL-2", "SELL", "SATB1", "BACH2", "S100A4", "S100A6",
    "ITGB1", "IKZF2", "FOXP3", "CD127lo", "CD25hi", "GITR", "Helios"
)
genes[which(!genes %in% rownames(treg))]

# Gene names that had to be fixed:
fixnames <- c(
    "CD95" = "FAS",
    "CD45RA" = "PTPRC",
    "L-selectin" = "SELL",
    "CTLA-4" = "CTLA4",
    "Ki67" = "MKI67",
    "BCL-2" = "BCL2",
    "CD127lo" = "IL7R",
    "CD25hi" = "IL2RA",
    "GITR" = "TNFRSF18",
    "Helios" = "IKZF2"
)
names(genes) <- genes
genes[names(fixnames)] <- fixnames
genes


csmarkers_plot <- patchwork::wrap_plots(
    show_tsneplot(treg),
    FeaturePlot(treg, features = csgenes, reduction = "tsne")
)
ggsave(fs::path(figdir, "clusters_raw.png"), show_tsneplot(treg), width = 5, height = 5)
ggsave(fs::path(figdir, "csmarkers.png"), FeaturePlot(treg, features = genes, reduction = "tsne"), width = 20, height = 25)
