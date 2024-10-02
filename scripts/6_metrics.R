# Scripts to provide summary statistics and metrics for the paper
source("code/metrics.R")
library(dplyr)

# How many patients in each group for TREG?
seu <- qs::qread(fs::path("gdm/results_v3/", "gdm_subtyped.qs"))$treg@meta.data
unique(seu$cluster)
seu %>%
    dplyr::select(multi_q, cond_cln) %>%
    unique() %>%
    dplyr::count(cond_cln)
# CONTROL: 13 ; GDM: 10 #
# Note in metadata, GDM is 9 as one sample is missing blood pressure.

# Mean and SD of cells per patient?
seu %>%
    dplyr::group_by(multi_q, cond_cln) %>%
    dplyr::transmute(ncells = n()) %>%
    unique() %>%
    ungroup() %>%
    dplyr::group_by(cond_cln) %>%
    dplyr::transmute(meanpp = mean(ncells), sdpp = sd(ncells)) %>%
    unique()

# Total number of cells per group?
seu %>%
    dplyr::select(cond_cln) %>%
    dplyr::count(cond_cln)

# Total cells?
dim(seu)

# What are the top-3 subtyping marker genes for each group?
ma <- readr::read_csv(fs::path("gdm/results_v2/", "top_subtyping_de_markers.csv")) %>%
    dplyr::select(ident, label, gene, cpva = CONTROL_p_val_adj, cal = CONTROL_avg_log2FC)

show_top_marks <- function(x, i, l = "treg", n = 3) {
    x %>%
        dplyr::filter(ident == i & label == l) %>%
        dplyr::top_n(n, abs(cal)) %>%
        dplyr::arrange(-abs(cal)) %>%
        dplyr::mutate(gsign = paste0(gene, ifelse(sign(cal) == 1, "+", "-"))) %>%
        pull(gsign) %>%
        paste0(collapse = ", ")
}
for (e in seq(0, 11)) {
    cat(paste0(show_top_marks(ma, e), "\n"))
}
for (e in seq(0, 10)) {
    cat(paste0(show_top_marks(ma, e, "cd4", 4), "\n"))
}


# View CCR7 expression
gdm <- qs::qread(fs::path(outdir, "gdm_seurat_sct_markers.qs"))
FeaturePlot(gdm$treg, "CCR7", reduction = "tsne")

# What are the propeller abundance estimates and p-values?
pp <- qs::qread(fs::path("gdm/results_v3/", "propeller.qs"))
tp <- pp[[1]]
cp <- pp[[2]]

# Patient cell proportions
abfigobs <- qs::qread("gdm/results_v3/ab_test_fig.qs")
xcp <- abfigobs$treg_ab[[4]]$data
dplyr::filter(xcp, cluster == "MALAT1+") %>% dplyr::arrange(-clust_prop)
xcd <- abtreg[[3]]$data
abfigobs$cd4_ab[[3]]$data

## Metrics on DEGs per cluster
dedata <- qs::qread(fs::path(outdir, "treg_fig3_obs.qs"))
dedata[[1]]$data %>%
    dplyr::group_by(label, cluster) %>%
    dplyr::summarise(sigc = sum(is_sig, na.rm = T), ng = n())

# Mean sig DE genes per cluster
dedata[[1]]$data %>%
    dplyr::group_by(label, cluster) %>%
    dplyr::summarise(sigc = sum(is_sig, na.rm = T), ng = n()) %>%
    dplyr::mutate(isnai = stringr::str_detect(cluster, "Naive")) %>%
    dplyr::group_by(isnai) %>%
    dplyr::mutate(msc = mean(sigc))

# Generate markerlist as table based on marker source
mkl <- qs::qread(fs::path(outdir, "markers.qs"))
markdf <- rbind(
    data.frame(
        gene = mkl$markerlist$clustmarks,
        source = "Treg subset marker gene"
    ),
    data.frame(
        gene = mkl$markerlist$de_marks,
        source = "Treg GDM differentially expressed gene"
    )
)
readr::write_csv(markdf, "gdm/data/marker_gene_stable.csv")

## GDM gene AUC
fca <- fs::path(outdir, "gdm_covar_figcache.qs")
fca <- qs::qread(fca)
bulkp <- fca[[4]][[1]]$data
pseup <- fca[[4]][[2]]$data

pseup %>% dplyr::filter(model_name %in% c("MT-CO3", "RPS18", "ANXA2"))
bulkp %>% dplyr::filter(model_name %in% c("HLA-DRA"))

bulkp %>%
    pull(model_name) %>%
    unique() %>%
    length() # N genes returned

top_markers <- readr::read_csv(fs::path(outdir, "top_subtyping_de_markers.csv"))

top_markers %>%
    dplyr::filter(gene == "TXNIP") %>%
    dplyr::select(-contains("CASE"))

bulkp %>% dplyr::filter(all_y_auc >= 0.65)
