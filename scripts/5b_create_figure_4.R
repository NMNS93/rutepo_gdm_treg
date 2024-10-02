# Create Figure 4 plots

# Setup
source("code/plots.R")
source("code/glm_meta_model.R")
library(data.table)
library(ggplot2)
library(ggupset)
outdir <- "gdm/results_v3/"
figdir <- fs::dir_create("gdm/figures_v3/covariates")
maindir <- fs::dir_create("gdm/figures_v3/main")

# A: Metadata
meta <- fread(fs::path(outdir, "meta_data.csv"))
p_meta <- plot_gdm_meta(meta)
odds_data <- data.table::fread(fs::path(outdir, "meta_odds.csv"))
p_odds <- plot_odds_rat(odds_data)
p_metas <- patchwork::wrap_plots(p_meta$indiv, nrow = 2)
# ggsave(fs::path(figdir, "metadata_all.png"), p_metas, width=8, height=8)

# B: Top single markers from all datasets
stirms_smod <- qs::qread(fs::path(outdir, "stirms_smod.qs"))
wang_smod <- qs::qread(fs::path(outdir, "wang_smod.qs"))
cd4_smod <- qs::qread(fs::path(outdir, "cd4_smod.qs"))
treg_smod <- qs::qread(fs::path(outdir, "treg_smod.qs"))


# Show markers
smod_all_tops <- qs::qread(fs::path(outdir, "smod_all_tops.qs"))

roc_data_af <- smod_simple_format(
    smod_all_tops[[2]], smod_all_tops[[1]], smod_all_tops[[3]], smod_all_tops[[4]]
)

roc_data_ss2 <- smod_to_roc_ss(
    smod_all_tops[[2]], smod_all_tops[[1]], smod_all_tops[[3]], smod_all_tops[[4]]
)

prd_f1_plots <- plot_roc_f1_single(roc_data_af)

# Export CD4 markers for supplementary table
c4tops <- smod_all_tops[[3]] %>%
    dplyr::select(data_name, model_name, all_y_auc) %>%
    unique() %>%
    dplyr::arrange(-all_y_auc) %>%
    dplyr::mutate(gene = stringr::str_remove(model_name, "status;")) %>%
    head(10)
data.table::fwrite(c4tops, "gdm/results_v3/c4_top_auc.csv")

# Show select genes that are in all datasets
roc_data_ss2 <- qs::qread(fs::path(outdir, "roc_data_ss2_cache.qs"))
rdsel <- roc_data_ss2 %>%
    dplyr::select(dataset, model_name) %>%
    unique() %>%
    dplyr::mutate(value = 1) %>%
    tidyr::pivot_wider(id_cols = "model_name", names_from = "dataset", values_from = "value")

drop_mods <- rdsel$model_name[(rdsel %>% is.na() %>% rowSums()) > 0]

top_rgenes_bulk <- get_top_rgenes_bulk(roc_data_ss2) %>% dplyr::filter(!model_name %in% drop_mods)
prd_roc_plot <- plot_roc_data_single(roc_data_ss2,
    rgenes = top_rgenes_bulk$model_name[1:6], ncol = 2
)


ggsave(fs::path(figdir, "roc_single.png"), prd_roc_plot, width = 4, height = 8)
ggsave(fs::path(figdir, "auc_f1_bulk.png"), prd_f1_plots[[1]], width = 8, height = 8)
ggsave(fs::path(figdir, "auc_f1_pseu.png"), prd_f1_plots[[2]], width = 8, height = 8)

# CD4/Treg top 5
## Load Treg top 10
tr5 <- qs::qread(fs::path(outdir, "treg_tops5.qs"))$tops
cup <- tr5 %>%
    dplyr::filter(all_y_f1 > .5) %>%
    dplyr::top_n(10, all_y_auc)
dcup <- tops5_to_upset(cup)

## Add CD4 scores at top 10
cd5all <- qs::qread(fs::path(outdir, "cd4_mod5.qs")) %>% # Long-running
    dplyr::filter(model_name %in% dcup$model_name) %>%
    dplyr::select(model_name, all_cd4_auc = all_y_auc) %>%
    unique()
dcupa <- dcup %>%
    dplyr::left_join(cd5all, by = c("model_name")) %>%
    tidyr::pivot_longer(c(all_y_auc, all_cd4_auc))

## Gnenerate plot
p_cd4_upset <- ggplot(dcupa, aes(x = genes, y = value)) +
    geom_col(
        aes(fill = name),
        position = position_dodge()
    ) +
    scale_x_upset() +
    theme_bw() +
    labs(x = "", y = "CD4+ T pseudobulk\nAUC (LOO-CV)") +
    scale_y_continuous(n.breaks = 10)
p_cd4_upset
ggsave(fs::path(figdir, "cd4_pseu_upset.png"), p_cd4_upset, width = 5, height = 8)
#

# RFE boxplots::::
rfe_res <- qs::qread(fs::path(outdir, "rfe_singles.qs"))
wang_rfe <- rfe_res[[1]]
stirm_rfe <- rfe_res[[2]]

wang_rfe_top <- wang_rfe[[1]]$optVariables
stirm_rfe_top <- stirm_rfe[[1]]$optVariables
wangl <- load_wang_bulk("gdm/data/GSE154414_Expression_Gene.csv")
stirm <- load_stirm_bulk("data/mrna/GSE92772_expression_mrna.tsv", "data/mrna/GSE92772_series_matrix.txt")

p_wang_rfe <- ggplot(
    dplyr::filter(wangl, gene %in% wang_rfe_top),
    aes(x = gene, y = log10(expression + 1), fill = status_c)
) +
    geom_boxplot(outlier.shape = NA) +
    labs(fill = "") +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "Wang et al. 2021.\nTop CD4+ Markers by RFE") +
    geom_jitter(position = position_jitterdodge(dodge.width = .75, jitter.width = 0.1, jitter.height = 0.1))

p_stirm_rfe <- ggplot(
    dplyr::filter(stirm, gene %in% stirm_rfe_top),
    aes(x = gene, y = expr, fill = status)
) +
    geom_boxplot(outlier.shape = NA) +
    labs(fill = "") +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "stirm et al. 2018.\nTop CD4+ Markers by RFE") +
    geom_jitter(position = position_jitterdodge(dodge.width = .75, jitter.width = 0.1, jitter.height = 0.1))

p_rfe <- patchwork::wrap_plots(p_wang_rfe, p_stirm_rfe, nrow = 1)
ggsave(fs::path(figdir, "plot_rfe_tops.png"), p_rfe, width = 10, height = 10)

# ## Figure 8 : DimReg on RFE ----

wdimreg <- dplyr::filter(wangl, gene %in% wang_rfe_top) %>% tidyr::pivot_wider(names_from = "gene", values_from = "expression")

p_wang_me <- ggplot(wdimreg, aes(x = RPL27A, y = TXNIP)) +
    geom_point(aes(fill = status_c), shape = 21, size = 7) +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Wang et al., 2021.\nTop Maker Expression")

wang_wide <- qs::qread(fs::path(outdir, "wang_wide.qs"))
l_p_wang <- plt_umap_dimreg(wang_wide, wang_rfe_top, 2)
p_wang_umap <- l_p_wang[[2]]
p_wang_umap

stw <- qs::qread(fs::path(outdir, "strims_wide.qs"))
l_p_stirm <- plt_umap_dimreg(stw, stirm_rfe_top, 5)
p_stirm_umap <- l_p_stirm[[2]]
p_stirm_umap
p_stirm_umap <- p_stirm_umap + labs(title = "Stirm et al., 2018.\nTop CD4+ Marker UMAP") + theme_bw() + theme(legend.position = "bottom")
p_dimreg <- patchwork::wrap_plots(p_wang_umap, p_stirm_umap, nrow = 1)

# Panel B - Pseudobulk Expression -----
cd4data <- qread(here::here("gdm/results_v3/cd4_model_dataset.qs"))[[2]]
c4plot <- ggplot(cd4data, aes(x = RPS18, y = `MT-CO3`)) +
    geom_point(aes(color = factor(status))) +
    theme_bw() +
    scale_color_manual(values = as.character(gdm_pal[c("orange", "black")]), labels = c("GDM", "Control")) +
    labs(color = "") +
    theme(legend.position = "none")
ggsave(fs::path(figdir, "c4_top2genes.png"), c4plot)

# Panel C - PCA on Wang and STIRM -----
mkl <- qs::qread(fs::path(outdir, "markers.qs"))$markers

# Wang
wangmark <- load_wang_bulk("gdm/data/GSE154414_Expression_Gene.csv") %>% filter(gene %in% mkl)
wamat <- tidyr::pivot_wider()

# Stirm
stmark <- load_stirm_bulk("data/mrna/GSE92772_expression_mrna.tsv", "data/mrna/GSE92772_series_matrix.txt") %>% filter(gene %in% mkl)
stmat <- tidyr::pivot_wider(stmark, id_cols = c(id, status), names_from = "gene", values_from = "expr")
stpc <- prcomp(stmat[, -1:-2], scale = F, center = T)
stpca <- ggplot(
    data.frame(pc1 = stpc$x[, 1], pc2 = stpc$x[, 2], status = stmat$status),
    aes(x = pc1, y = pc2, color = status)
) +
    geom_point() +
    theme_bw()
ggsave(fs::path(figdir, "stirm_pca_all.png"), stpca) # Poor separation

#  (Described in text only) - CD4 AUC scores ----
cd4_smod <- qs::qread(fs::path(outdir, "cd4_smod.qs"))
cmod_ranks <- cd4_smod %>%
    roc_mini() %>%
    dplyr::select(dataset, model_name, all_y_auc) %>%
    unique() %>%
    dplyr::group_by(model_name) %>%
    dplyr::summarise(mauc = median(all_y_auc)) %>%
    dplyr::filter(mauc >= 0.6) %>%
    dplyr::mutate(model_name = stringr::str_remove(model_name, "status;")) %>%
    dplyr::arrange(-mauc)
cmod_ranks$model_name <- factor(cmod_ranks$model_name, levels = cmod_ranks$model_name)
p_cmod <- ggplot(cmod_ranks, aes(y = model_name, x = mauc)) +
    geom_col() +
    theme_bw()
ggsave(fs::path(figdir, "cd4_smod_auc.png"), p_cmod)

#  AUC and F1 scores ----
# Single marker scores from both datasets
stirms_smod <- qs::qread(fs::path(outdir, "stirms_smod.qs"))
wang_smod <- qs::qread(fs::path(outdir, "wang_smod.qs"))
afscores <- rbind(
    dplyr::select(stirms_smod, model_name, data_name, all_y_auc, all_y_f1),
    dplyr::select(wang_smod, model_name, data_name, all_y_auc, all_y_f1)
) %>%
    unique() %>%
    dplyr::mutate(model_name = str_remove(model_name, "status;"))
afscores$datan2 <- ifelse(afscores$data_name == "stirms_wide", "Stirm et al.", "Wang et al.")
### cache for supplementary table
data.table::fwrite(dplyr::filter(afscores, all_y_auc >= 0.7), "/genomics/projects/B061-PanicosShangaris/nana/gdm/results_v3/stirm_wang_top_supplement.csv")
txtdf <- dplyr::filter(afscores, all_y_auc >= 0.7)
afplot <- ggplot(afscores, aes(x = all_y_f1, y = all_y_auc)) +
    geom_hline(yintercept = 0.71, linetype = "dashed") +
    geom_jitter(aes(color = datan2, shape = datan2), alpha = 0.7, size = ifelse(afscores$data_name == "stirms_wide", 4, 3.5), width = 0.005, height = 0.005) +
    theme_bw() +
    geom_label_repel(data = txtdf, aes(label = model_name, color = datan2), max.overlaps = 10, size = 3, nudge_x = -.24, nudge_y = .06, box.padding = 0.1, show_guide = F) +
    labs(x = "F1 Score", y = "AUC") +
    scale_color_manual(values = as.character(gdm_pal[c("lightpurple", "teal")])) +
    labs(color = "Study", shape = "Study") +
    scale_y_continuous(limits = c(0.5, NA)) +
    theme(legend.position = c(0.1, .15))
ggsave(fs::path(figdir, "wang_stirms_auc_f1.png"), afplot)

# Raw expression values -----
wange <- load_wang_bulk("gdm/data/GSE154414_Expression_Gene.csv") %>%
    filter(gene %in% (txtdf %>% dplyr::filter(data_name == "wang_wide") %>% dplyr::top_n(5, all_y_auc) %>% pull(model_name))) %>%
    dplyr::select(id, gene, expr = expression, status = status_c) %>%
    mutate(Study = "Wang et al.")
stirme <- load_stirm_bulk("data/mrna/GSE92772_expression_mrna.tsv", "data/mrna/GSE92772_series_matrix.txt") %>%
    filter(gene %in% (txtdf %>% dplyr::filter(data_name == "stirms_wide") %>% dplyr::top_n(5, all_y_auc) %>% pull(model_name))) %>%
    dplyr::mutate(status = ifelse(status == "GDM", "GDM", "CONTROL"), Study = "Stirm et al.")
rawe <- rbind(wange, stirme) %>%
    mutate(`log10(Expression+1)` = log10(expr + 1)) %>%
    dplyr::filter(gene != "MT-ATP8")
p_expr2 <- ggplot(rawe, aes(y = gene, x = `log10(Expression+1)`, color = status)) +
    geom_point(aes(group = status), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
    scale_color_manual(values = as.character(gdm_pal[c("darkgrey", "orange")])) +
    stat_summary(
        aes(group = status),
        fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", width = 0.5, lwd = 0.2,
        # add this bit here to your stat_summary function
        position = position_dodge(width = 0.8)
    ) +
    facet_wrap(~Study, nrow = 1, scales = "free") +
    theme_bw() +
    labs(y = "") +
    theme(legend.position = "bottom", strip.background = element_blank())
ggsave(fs::path(figdir, "gene_expr_wang_stirm.png"), p_expr2)

# TXNIP expression in single cell datasets -----
gdm <- qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 <- gdm$cd4
treg <- gdm$treg
rm(gdm)
# treg_txnip = FeaturePlot(treg, features = c("TXNIP"), split.by = "control")
Idents(treg) <- "cluster"
treg_txnip <- VlnPlot(treg, "TXNIP", split.by = "control", log = F, add.noise = T, slot = "scale.data", pt.size = 0.00005, cols = c("orange", "lightgrey")) + labs(title = "Treg | TXNIP")
Idents(cd4) <- "cluster"
cd4_txnip <- VlnPlot(cd4, "TXNIP", split.by = "control", log = F, add.noise = T, slot = "scale.data", pt.size = 0.00005, cols = c("orange", "lightgrey")) + labs(title = "CD4+ | TXNIP")
cpatch <- patchwork::wrap_plots(treg_txnip, cd4_txnip, ncol = 1) +
    plot_annotation(tag_levels = "A")
ggsave(filename = fs::path(figdir, "TXNIP.png"), plot = cpatch, width = 8, height = 8)


# Figures ------------
#
# Cache figure/plot data in case require replotting
figcache <- list(
    # p_metas,
    p_odds, # USE
    p_cd4_upset,
    prd_roc_plot,
    prd_f1_plots,
    p_rfe,
    p_dimreg,
    c4plot, # USE
    afplot, # USE
    p_expr2 # USE
)
qs::qsave(figcache, fs::path(outdir, "gdm_covar_figcache_v2.qs"))

# # Create Figure panel ----

# Wrap_elements must be used to be able to label with patchwork
# and to ensure the upset plot x axis does not override
p_panel2 <- wrap_elements(full = p_odds) +
    wrap_elements(full = c4plot) +
    wrap_elements(full = afplot) +
    wrap_elements(full = p_expr2) +
    plot_layout(design = "
              ABDDD
              CCDDD
              CCDDD
              ") +
    plot_annotation(tag_levels = c("A"))

# Note: Cowplot required because upsetplot errors
ggsave(fs::path(maindir, "Figure4_v2.png"), p_panel2, width = 12, height = 8)
qs::qsave(p_panel2, fs::path(outdir, "gdm_fig4_obs.qs"))
