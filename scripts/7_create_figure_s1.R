# Compose Supplementary Figure 1: CD4+ T cells
source("code/plots.R")
source("code/subtyping_de.R")
library(qs)

outdir <- here::here("gdm/results_v3")
figdir <- fs::dir_create("gdm/figures_v3/main")
gdm <- qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 <- gdm$cd4
treg <- gdm$treg

## Recapitulate all_cd4_fig1
# Plot TSNE for CD4+ T
Idents(cd4) <- "cluster"
cd4_tsne <- show_tsneplot(cd4, pal = cd4_dimreg_pal) + labs(title = "CD4+ T Cells")

## Show highlight features
cd4_highlight <- c("CCR7", "IL2RA", "CCR6", "CXCR3", "S100A4", "FOXP3")
p_cd4high <- FeaturePlot(cd4, features = cd4_highlight, cols = c("lightgray", gdm_pal[["orange"]]), reduction = "tsne") + NoLegend()

# Show heatmap for CD4
top_markers <- readr::read_csv(fs::path("gdm/results_v2/", "top_subtyping_de_markers.csv"))
markers <- de_top_n_genes(dplyr::filter(top_markers, label == "cd4"), 2)
markers <- c(markers, "FOXP3", cd4_highlight)
cd4_avc <- qs::qread(fs::path(outdir, "gdm_average_expr.qs"))$cd4_avc
cd4_hm <- subtype_heatmap(cd4_avc, markers, is.treg = F)

# Show humanin+ boxplots for CD4 and Treg
library(Seurat)
library(ggplot2)
library(dplyr)

get_data <- function(seurat_obj, sample_name, feature) {
    data.frame(
        sample = sample_name,
        count = seurat_obj@assays$RNA@data[feature, ],
        cluster = seurat_obj@meta.data$cluster,
        feature = feature
    )
}
fdata <- data.table::rbindlist(list(
    get_data(cd4, "CD4", "MTRNR2L12"),
    get_data(cd4, "CD4", "MTRNR2L8"),
    get_data(cd4, "CD4", "MALAT1"),
    get_data(treg, "Treg", "MTRNR2L12"),
    get_data(treg, "Treg", "MTRNR2L8"),
    get_data(treg, "Treg", "MALAT1")
))
fdata$cluster <- dplyr::case_when(
    fdata$cluster == "Humanin+" ~ "Humanin+",
    fdata$cluster == "MALAT1+" ~ "MALAT1+",
    TRUE ~ "Other"
)

p_humanin <- ggplot(fdata, aes(x = cluster, y = count, fill = sample)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    labs(
        fill = "Dataset",
        y = "Feature count",
        x = "Cluster"
    ) +
    scale_fill_manual(values = as.character(gdm_pal[c("teal", "orange")])) +
    facet_wrap(~feature) +
    theme(axis.text.x = element_text(angle = 90))

## Plot differential abundance
pp <- qs::qread(fs::path(outdir, "propeller.qs"))
tpropeller <- pp[[1]]
cpropeller <- pp[[2]]
p_ab_test_cd4 <- plot_propeller_res(cpropeller %>% filter(cluster == "Memory-1"), "CD4", gdm_pal[["pink"]])
p_ab_test_cd4

## Compose plot
qs::qsave(list(cd4_tsne, p_cd4high, cd4_hm, p_humanin), fs::path(outdir, "supp_fig.qs"))
p_supp_fig <- patchwork::wrap_plots(
    cd4_tsne, p_cd4high, cd4_hm, p_humanin
) + plot_layout(design = "
AAABBB
AAABBB
CCCDDD
CCCDDD
") + plot_annotation(tag_levels = list(c("A", "B", "", "", "", "", "", "C", "D")))
ggsave(fs::path(figdir, "supp_fig1.png"), p_supp_fig, width = 14, height = 12, dpi = 300)


## Supplementary Table. ----
## Shows the marker genes used for GDM analysis and their origin.

## Get a unique list of marker genes
mkl <- qs::qread(fs::path("gdm/results_v2", "markers.qs"))
markers <- mkl$markers
marklist <- mkl$markerlist

## For marker gene sources, setup a look up table for their supplement name
codes <- c("rdeg_treg", "rdeg_cd4", "nact_treg", "hub_string", "humanin")
stopifnot(all(codes %in% names(marklist)))
code_rename <- c("Treg DEG", "CD4+ DEG", "Treg Naive-Act-2 DEG", "Treg Naive-Act-2 Hub Gene", "Treg Humanin+ DEG")
names(code_rename) <- codes
## Write a function to get marker gene sources per gene
gmgspg <- function(x) {
    cond <- sapply(marklist, function(i) x %in% i)
    mcodes <- codes[which(cond)]
    if (length(mcodes) == 0) {
        stop("Error: all markers should have codes.")
    }

    mcln <- as.character(
        code_rename[mcodes]
    )
    mcln <- paste0(mcln, collapse = ";")
    return(mcln)
}


## Create marker gene results table
mgt <- data.frame(
    "Gene" = markers,
    "Criteria" = sapply(markers, function(e) gmgspg(e))
)
data.table::fwrite(mgt, "gdm/results_v2/marker_genes.csv")
