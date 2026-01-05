# Perform differential expression testing between clusters to determine which genes identify each cluster.

# Library ----
source("code/helpers.R")
source("code/plots.R")
source("code/subtyping_de.R")
library(dplyr)

# Variables ----
outdir <- fs::dir_create("gdm/results/")
figdir <- fs::dir_create("gdm/figures/subtyping")
gdm <- qs::qread(fs::path(outdir, "gdm_seurat_reclust_filter.qs"))

# Main ----

# Setup Seurat datasets
treg <- gdm$treg
cd4 <- gdm$cd4
rm(gdm)
DefaultAssay(treg) <- "SCT"
DefaultAssay(cd4) <- "SCT"

# Prepare samples for DE testing
treg <- PrepSCTFindMarkers(treg)
cd4 <- PrepSCTFindMarkers(cd4)

# Run differential expression
treg_markers <- subtyping_de(treg, "treg")
cd4_markers <- subtyping_de(cd4, "cd4")
subtyping_de_markers <- rbind(treg_markers, cd4_markers)
readr::write_csv(subtyping_de_markers, fs::path(outdir, "subtyping_de_markers.csv"))

# Select top markers
top_markers <- filter_for_top_markers(subtyping_de_markers)
readr::write_csv(top_markers, fs::path(outdir, "top_subtyping_de_markers.csv"))
qs::qsave(list("treg" = treg, "cd4" = cd4), fs::path(outdir, "gdm_seurat_sct_markers.qs"))

# Exploratory visualisation (Treg)
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

# Exploratory visualisation (CD4)
de_cd4 <- readr::read_csv(fs::path(outdir, "top_subtyping_de_markers.csv")) %>% dplyr::filter(label == "cd4")
p_de_cd4_top3 <- ShowSCTClusterDeTop(cd4, de_cd4)
cd4_de_idents <- seq_along(unique(de_cd4$ident))
for (i in cd4_de_idents) {
  ggsave(fs::path(figdir, paste0("p_cd4_subtype_top3_", i, ".png")), p_de_cd4_top3[[i]], width = 15, height = 6)
}
