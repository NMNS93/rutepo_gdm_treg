library(Seurat)
library(scRepertoire)
library(dplyr)
library(ggplot2)
library(data.table)
library(parallel)
source("code/preprocessing_qc.R")

# Variables
outdir <- fs::dir_create(here::here("data/seurat/220802_Shangaris"))
figdir <- fs::dir_create(here::here("figures/220802_seurat"))
meta_csv <- here::here("data/metadata/220721_metadata_cln.csv")
tenx <- qs::qread(fs::path(outdir, "seurat_tenx_clustered.qs"))
DefaultAssay(tenx) <- "RNA"
tenx <- AddPanicosMeta(tenx, meta_csv, cols_regex = "developed|condition|control")
cl <- qs::qread(fs::path(outdir, "screp_combined_tcr.qs")) # clonotypes

# Reorder
# As the data is plotted in the order of the list, here we create a metadata table to
# alter the plots downstream
ord <- data.table(samples = names(cl))
mng <- tidyr::separate(ord, samples, c("patient", "condition"), sep = "_", remove = F, extra = "merge") %>%
    tidyr::separate(patient, into = c("run", "hashtag"), sep = -5, remove = F) %>%
    dplyr::mutate(run = str_remove(run, "-$")) %>%
    tidyr::separate(run, into = c("base", "ix", "ctype"), sep = "-", remove = F) %>%
    dplyr::arrange(ctype, condition, base, hashtag)
cl <- cl[mng$samples]

# Data
# Fist, we quantify the number of unique clonotypes across samples
# This function is a bar plot over the CT_strict field in each dataset
## Unique clonotypes may represent batch effect. Unclear how scale=T improves result.
v <- quantContig(cl, cloneCall = "strict", scale = F) +
    coord_flip() + theme(legend.position = "none") + scale_fill_manual(values = rep("red", 58))
p_uclon <- (v$data) %>%
    left_join(dplyr::rename(mng, "values" = "samples")) %>%
    ggplot(aes(y = values, x = total, fill = condition)) +
    geom_col() +
    theme_bw() +
    facet_wrap(~ctype, scales = "free_y") +
    labs(x = "unique_clonotypes") +
    theme(text = element_text(size = 8))
ggsave(here::here("figures/unique_clonotypes.png"), p_uclon, width = 8, height = 5)

# Indeed, contig abundance demonstrates that most clonotypes are present in few.
p_abu <- abundanceContig(cl) + theme(legend.position = "none") + scale_color_manual(values = rep("black", length(cl)))
ggsave(here::here("figures/abundance_clonotypes.png"), p_abu, width = 8, height = 5)

# For patients in each condition, plot the proportion of the most abundant clonotype
conds <- c("CONTROL", "GDM", "PE")
plt_cond <- function(cl, x, ctype = "Treg") {
    pats <- names(cl) %>% purrr::keep(~ stringr::str_detect(.x, paste0(ctype, ".*", x)))
    compareClonotypes(cl,
        samples = pats,
        numbers = 1,
        cloneCall = "aa",
        graph = "alluvial"
    ) + theme(legend.position = "none", text = element_text(size = 5)) +
        coord_flip() +
        labs(title = paste0(ctype, "_", x))
}

p_conds_tregs <- lapply(conds, function(x) {
    plt_cond(cl, x)
})
p_conds_cd4 <- lapply(conds, function(x) {
    plt_cond(cl, x, ctype = "CD4")
})

ggsave(here::here("figures/p_conds_tregs.png"), patchwork::wrap_plots(p_conds_tregs, nrow = 1), width = 8, height = 5)
ggsave(here::here("figures/p_conds_cd4.png"), patchwork::wrap_plots(p_conds_cd4, nrow = 1), width = 8, height = 5)

# Visualise gene usage
plt_genesViz <- function(cl, x, ctype) {
    pats <- names(cl) %>%
        purrr::keep(~ stringr::str_detect(.x, paste0(ctype, ".*", x))) %>%
        purrr::discard(~ stringr::str_detect(.x, "L226-A-Treg-AHH05_PE"))
    vizGenes(cl[pats],
        gene = "V",
        chain = "TRB",
        y.axis = "J",
        plot = "heatmap",
        scale = TRUE,
        order = "gene"
    ) + labs(title = paste0(ctype, "_", x))
}

p_gV_tregs <- lapply(conds, function(x) {
    plt_genesViz(cl, x, "Treg")
})
p_gV_cd4 <- lapply(conds, function(x) {
    plt_genesViz(cl, x, ctype = "CD4")
})

ggsave(here::here("figures/p_gv_tregs.png"), patchwork::wrap_plots(p_gV_tregs, nrow = 1), width = 15, height = 5)
ggsave(here::here("figures/p_gv_cd4.png"), patchwork::wrap_plots(p_gV_cd4, nrow = 1), width = 15, height = 5)

# Filter: Remove sample without clonotype info
mng <- dplyr::filter(mng, samples != "L226-A-Treg-AHH05_PE")
cl <- cl[mng$samples]

plt_homeo <- function(cl, x) {
    pats <- names(cl) %>%
        purrr::keep(~ stringr::str_detect(.x, paste0("[0-9]+_", x))) %>%
        purrr::discard(~ stringr::str_detect(.x, "L226-A-Treg-AHH05_PE"))
    clonalHomeostasis(cl[pats],
        cloneCall = "gene",
        cloneTypes = c(
            Rare = 1e-04,
            Small = 0.001,
            Medium = 0.01,
            Large = 0.1,
            Hyperexpanded = 1
        )
    ) +
        labs(title = x) + coord_flip()
}
p_homeo <- patchwork::wrap_plots(lapply(conds, plt_homeo, cl = cl), ncol = 1)
ggsave(here::here("figures/vdj_homestasis.png"), p_homeo, width = 8, height = 8)

plt_coverlap <- function(cl, x) {
    pats <- names(cl) %>% purrr::keep(~ stringr::str_detect(.x, paste0("[0-9]+_", x, "$")))
    clonalOverlap(cl[pats],
        cloneCall = "gene+nt",
        method = "jaccard"
    ) +
        theme(axis.text.x = element_text(angle = 90), text = element_text(size = 8)) +
        labs(title = x)
}
p_coverlap <- lapply(conds, plt_coverlap, cl = cl)
for (i in seq_along(p_coverlap)) {
    ggsave(here::here(paste0("figures/vdj_clonalOverlaps", i, ".png")), p_coverlap[[i]], width = 8, height = 5)
}

# Diversity
cl2 <- lapply(
    cl,
    function(x) {
        dplyr::left_join(x, dplyr::rename(mng, "sample" = "patient")) %>%
            dplyr::filter(condition %in% c("PE", "GDM", "CONTROL"))
    }
)

cdivers_Treg <- clonalDiversity(cl2[names(cl2) %>% purrr::keep(~ stringr::str_detect(.x, "Treg"))],
    cloneCall = "gene",
    group.by = "samples",
    x.axis = "condition",
    n.boots = 100
) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position = "none") +
    labs(title = "Treg")
cdivers_CD4 <- clonalDiversity(cl2[names(cl2) %>% purrr::keep(~ stringr::str_detect(.x, "CD4"))],
    cloneCall = "gene",
    group.by = "samples",
    x.axis = "condition",
    n.boots = 100
) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position = "none") +
    labs(title = "CD4")
ggsave(here::here(paste0("figures/clonalDiversity.png")),
    patchwork::wrap_plots(cdivers_Treg, cdivers_CD4, nrow = 1),
    width = 14, height = 5
)
