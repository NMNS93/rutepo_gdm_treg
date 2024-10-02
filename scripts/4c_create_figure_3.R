source("code/plots.R")

# Setup
library(tidyHeatmap)
library(tibble)
library(dplyr)
library(ggplot2)
figdir <- fs::dir_create("gdm/figures_v3/main/")
outdir <- "gdm/results_v3"
# Load DE objects
de_ <- qs::qread(fs::path(outdir, "we_treg_cln.qs"))
de_treg <- de_$wta
de_cd4 <- de_$wca
rm(de_)
de_treg <- dplyr::filter(de_treg, !(seurat_clusters %in% c(5, 6))) # Remove low cell clusters
# Genes sig DE in any cluster

# Load GSEA objects
tga <- qs::qread(fs::path(outdir, "gsea_treg.qs"))
tga <- dplyr::filter(tga, !(seurat_clusters %in% c(5, 6))) # Remove low cell clusters

# Plot DE dots
de_treg$is_sig <- de_treg$is_signif
# Set genes to those significant in eatleast one cluster
tgenes <- de_treg %>%
    dplyr::filter(is_signif == T) %>%
    pull(gene) %>%
    unique()
we_treg <- dplyr::filter(de_treg, gene %in% tgenes)
# Set levels of gene so that exclusive genes are ranked
glevels <- levels((we_treg %>% group_by(cluster) %>%
    mutate(gene = factor(gene, levels = gene[is_signif])) %>%
    ungroup())$gene)
we_treg$gene <- forcats::fct_relevel(we_treg$gene, glevels)
p_de <- de_dot_plot(
    we_treg,
    col.scale = c(gdm_pal[["ogreen"]], "lightgray", gdm_pal[["orange"]])
)

# Plot hallmark GSEA
tgap <- tibble(tga) %>%
    dplyr::mutate(
        sig_stars = dplyr::case_when(
            padj <= 0.001 ~ "***",
            padj <= 0.01 ~ "**",
            padj <= 0.05 ~ "*",
            padj <= 0.1 ~ "-",
            TRUE ~ ""
        ),
        n_edge = sapply(leadingEdge, length),
        pathway_cln =
            stringr::str_replace_all(
                stringr::str_remove(pathway, "HALLMARK_"),
                "_", " "
            )
    )

p_NES <- heatmap(
    tibble(tgap),
    cluster, pathway_cln, NES,
    scale = "row",
    palette_value = gdm_pal[c("ogreen", "white", "orange")],
    rect_gp = grid::gpar(col = "#161616", lwd = 0.5)
) %>%
    layer_text(.value = sig_stars, .size = 10) %>%
    wrap_heatmap(
        padding = grid::unit(c(0, -80, 0, -100), "points"), # Aesthetics w/ wrap heatmap b, l, t, r
        clip = FALSE
    )

qs::qsave(list(p_de, p_NES), fs::path(outdir, "treg_fig3_obs.qs"))

# Create fig 3.
treg_f3 <- patchwork::wrap_plots(
    p_de, p_NES,
    ncol = 1,
    heights = c(1, .75)
) + plot_annotation(
    tag_levels = list(c("A", "B"))
)

ggsave(fs::path(figdir, "Figure3.png"), treg_f3, width = 8, height = 12)
