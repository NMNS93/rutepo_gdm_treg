# Create main manuscript figures from existing data objects

# Global setup ----
library(ggplot2)
library(data.table)
source("code/plots.R")
outdir <- "gdm/results"
datadir <- "gdm/data"
figdir <- fs::dir_create("gdm/figures/main")
zdir <- fs::dir_create("gdm/export/") # Exports datasets underlying figures

# Figure 1 ----

## Library ----
source("code/plots.R")
source("code/subtyping_de.R")
library(data.table)
library(ggplot2)

## Variables ----
gdm <- qs::qread(fs::path(outdir, "gdm_subtyped.qs"))

## Main ----

# Load Seurat RNA datasets
cd4 <- gdm$cd4
treg <- gdm$treg
rm(gdm)

# Load GDM expression averages
saverages <- qs::qread(fs::path(outdir, "gdm_average_expr.qs"))
treg_avc <- saverages$treg_avc
cd4_avc <- saverages$cd4_avc
rm(saverages)

# Load DE marker lists
treg_top <- readr::read_csv(fs::path(outdir, "top_subtyping_de_markers.csv")) %>%
  dplyr::filter(label == "treg")
cmarkers <- readr::read_csv(fs::path(datadir, "paper_markers.csv"))

## Plots -----

# A
Idents(treg) <- "cluster"
treg_tsne <- show_tsneplot(treg, pal = treg_dimreg_pal) + labs(title = "Treg")

# B
treg_highlight <- cmarkers$show[!is.na(cmarkers$show)]
p_trhigh <- FeaturePlot(treg,
  features = treg_highlight,
  cols = c(alpha("lightgray", .15), gdm_pal[["orange"]]), reduction = "tsne", ncol = 3
) &
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank())
p_trhigh

# C
top_n <- de_top_n_genes(treg_top, 5)
top_papers <- cmarkers$genes
genelist <- unique(c(top_n, top_papers))

# Set cluster ordering for heatmap
dataset <- treg_avc@assays$SCT@scale.data
subtype_order <- treg_avc@meta.data %>%
  dplyr::select(cluster, clust_group) %>%
  unique() %>%
  dplyr::arrange(clust_group) %>%
  pull(cluster)
genelist_ord <- reorder_by_expr(genelist, treg_avc@assays$SCT@scale.data, (subtype_order), 10)

# Get list of marker genes for heatmap
treg_hm <- subtype_heatmap(treg_avc, genelist_ord, is.treg = T)

# Render figure
treg_f1 <- patchwork::wrap_plots(
  treg_tsne, p_trhigh, treg_hm + coord_flip(),
  nrow = 2,
  heights = c(1, 1, 0.5),
  design = "AAABBB
          AAABBB
          CCCCCC
          CCCCCC"
) + plot_annotation(
  tag_levels = list(c("A", "B", rep("", 8), "C"))
)

ggsave(fs::path(figdir, "Figure1.png"), treg_f1, width = 12, height = 10)

## Export ----

data.table::fwrite(
  treg_tsne$data,
  fs::path(zdir, "f1_treg_tsne.csv")
)

#' Helper. Reformat marker gene TSNE (B) for export
reformat_p_tr_high_for_export <- function(x) {
  x <- x$data
  gene <- names(x)[[4]]
  names(x)[[4]] <- "fill"
  x["gene"] <- gene
  return(x)
}
b_data <- data.table::rbindlist(
  lapply(p_trhigh, reformat_p_tr_high_for_export)
)
data.table::fwrite(
  b_data,
  fs::path(zdir, "f1_treg_tsne_genes.csv")
)

data.table::fwrite(
  treg_hm$data,
  fs::path(zdir, "f1_treg_heatmap.csv")
)

# Figure 2 -----

## Variables ----
abfigobs <- qs::qread("gdm/results/ab_test_fig.qs") # Abundance figure objects
ctf <- fread(fs::path("data", "cytof", "percentages_clusters_Cytof_forNana.csv")) # Cytof data

## Main ----

# Parse abundance figure objects
abtreg <- abfigobs$treg_ab

# Set treg pallette
f1pal_order <- c("orange", "darkred", "blue", "ogreen", "darkgrey", "purple", "pink")
f2pal <- gdm_dimreg_pal[match(f1pal_order, names_gdm_dimreg_pal)]

# Format figure
abtreg_2 <- abtreg[[2]] + scale_fill_manual(values = f2pal) +
  labs(fill = "scRNA\ncluster", title = "Treg\nscRNA") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45))

# Format cytof data
ctf <- ctf %>%
  tidyr::pivot_longer(cols = c("per_GDM", "per_HC")) %>%
  dplyr::mutate(group = ifelse(name == "per_GDM", "GDM", "CONTROL")) %>%
  dplyr::mutate(value = value / 100)

# Set Cytof palette
cfp_colors <- c(
  gdm_pal[["black"]],
  gdm_pal[["lightgrey"]],
  gdm_pal[["brown"]],
  colorspace::lighten(gdm_pal[["brown"]], 0.9),
  colorspace::lighten(gdm_pal[["brown"]], 0.7),
  colorspace::lighten(gdm_pal[["brown"]], 0.5),
  colorspace::lighten(gdm_pal[["brown"]], 0.3),
  colorspace::lighten(gdm_pal[["orange"]], 0),
  colorspace::lighten(gdm_pal[["orange"]], 0.5),
  colorspace::lighten(gdm_pal[["purple"]], 0.1),
  colorspace::lighten(gdm_pal[["purple"]], 0.3),
  colorspace::lighten(gdm_pal[["pink"]], 0),
  colorspace::lighten(gdm_pal[["ogreen"]], 0.01),
  colorspace::lighten(gdm_pal[["ogreen"]], 0.5)
)

# Plot Cytof data
cfp <- ggplot(ctf, aes(x = group, y = value, fill = cluster)) +
  geom_col(color = "black") +
  labs(y = "Proportion", title = "CD4+ T\nCytof", x = "", fill = "Cytof\ncluster") +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
  scale_fill_manual(values = cfp_colors)

# Select additional plots
patplot <- abtreg[[4]] + labs(x = "Patient", fill = "Disease\nStatus")
abplt <- abtreg[[3]] + labs(x = "Ratio\n(GDM:CONTROL)")

# Render figure 2
treg_f2 <- patchwork::wrap_plots(
  abtreg_2,
  patchwork::free(patplot, type = "label"),
  patchwork::free(abplt, type = "label"),
  cfp,
  nrow = 1,
  widths = c(0.2, 1, 0.4, 0.2)
) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") & theme(
  legend.position = "bottom",
  legend.direction = "horizontal", # Make legend horizontal
  legend.title = element_text(hjust = 0.5), # Center title (optional)
  legend.box = "vertical", # Stack title above legend
  plot.margin = margin(5, 5, 5, 5)
)
ggsave(fs::path(figdir, "Figure2_cytof.png"), treg_f2, width = 10, height = 8.4)

## Export ----

f2_a <- dplyr::left_join(
  (abtreg_2$data %>% dplyr::count(cond_cln, cluster)),
  (abtreg_2$data %>% dplyr::count(cond_cln) %>% dplyr::rename(total = n))
) %>%
  dplyr::mutate(
    proportion = n / total,
    label = "treg"
  ) %>%
  dplyr::rename(condition = cond_cln)
data.table::fwrite(f2_a, fs::path(zdir, "f2_treg_proportions.csv"))

data.table::fwrite(
  patplot$data %>% dplyr::rename(
    patient = multi_q,
    condition = cond_cln
  ),
  fs::path(zdir, "f2_treg_proportions_patient.csv")
)

data.table::fwrite(abplt$data, fs::path(zdir, "f2_propeller_results.csv"))

data.table::fwrite(cfp$data, fs::path(zdir, "f2_cytof_proportions.csv"))

# Figure 3 ----

## Library ----
source("code/plots.R")
library(tidyHeatmap)
library(tibble)
library(dplyr)
library(grid)
library(ggplot2)

## Variables ----
figdir <- fs::dir_create("gdm/figures/main/")
outdir <- "gdm/results"
de_ <- qs::qread(fs::path(outdir, "we_treg_cln.qs"))
de_treg <- de_$wta
de_cd4 <- de_$wca
rm(de_)
de_treg <- dplyr::filter(de_treg, !(seurat_clusters %in% c(5, 6))) # Remove low cell clusters

## Main ----

# Load GSEA objects
tga <- qs::qread(fs::path(outdir, "gsea_treg.qs"))
tga <- dplyr::filter(tga, !(seurat_clusters %in% c(5, 6))) # Remove low cell clusters

# Set genes to those significant in eatleast one cluster
de_treg$is_sig <- de_treg$is_signif
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
  column_title_gp = gpar(fontsize = 0),
  row_title_gp = gpar(fontsize = 0),
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

ggsave(fs::path(figdir, "Figure3_v2.png"), treg_f3, width = 8, height = 12)

## Export ----
f3_obs <- qs::qread(fs::path(outdir, "treg_fig3_obs.qs"))
data.table::fwrite(f3_obs[[1]]$data, fs::path(zdir, "f3_differential_expression.csv"))
data.table::fwrite(tgap, fs::path(zdir, "f3_gsea.csv"))

# Figure 4 ----

## Library ----
source("code/plots.R")
source("code/glm_meta_model.R")
library(data.table)
library(ggplot2)
library(ggupset)

## Variables ----
outdir <- "gdm/results/"
figdir <- fs::dir_create("gdm/figures/covariates")
maindir <- fs::dir_create("gdm/figures/main")

## Main ----

# A: Metadata
meta <- fread(fs::path(outdir, "meta_data.csv"))
p_meta <- plot_gdm_meta(meta)
odds_data <- data.table::fread(fs::path(outdir, "meta_odds.csv"))
p_metas <- patchwork::wrap_plots(p_meta$indiv, nrow = 2)
p_odds <- plot_odds_rat(odds_data) # Main plot

# B: Top single markers from all datasets
stirms_smod <- qs::qread(fs::path(outdir, "stirms_smod.qs"))
wang_smod <- qs::qread(fs::path(outdir, "wang_smod.qs"))
cd4_smod <- qs::qread(fs::path(outdir, "cd4_smod.qs"))
treg_smod <- qs::qread(fs::path(outdir, "treg_smod.qs"))
smod_all_tops <- qs::qread(fs::path(outdir, "smod_all_tops.qs"))
roc_data_af <- smod_simple_format(
  smod_all_tops[[2]], smod_all_tops[[1]], smod_all_tops[[3]], smod_all_tops[[4]]
)
roc_data_ss2 <- smod_to_roc_ss(
  smod_all_tops[[2]], smod_all_tops[[1]], smod_all_tops[[3]], smod_all_tops[[4]]
)

# Exploratory ROC F1 plot
prd_f1_plots <- plot_roc_f1_single(roc_data_af)

# Export CD4 markers for supplementary table
c4tops <- smod_all_tops[[3]] %>%
  dplyr::select(data_name, model_name, all_y_auc) %>%
  unique() %>%
  dplyr::arrange(-all_y_auc) %>%
  dplyr::mutate(gene = stringr::str_remove(model_name, "status;")) %>%
  head(10)
data.table::fwrite(c4tops, "gdm/results_v3/c4_top_auc.csv")

# Select genes that are in all datasets
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

# Exploratory intermediate figures
ggsave(fs::path(figdir, "roc_single.png"), prd_roc_plot, width = 4, height = 8)
ggsave(fs::path(figdir, "auc_f1_bulk.png"), prd_f1_plots[[1]], width = 8, height = 8)
ggsave(fs::path(figdir, "auc_f1_pseu.png"), prd_f1_plots[[2]], width = 8, height = 8)

# CD4/Treg top 5
tr5 <- qs::qread(fs::path(outdir, "treg_tops5.qs"))$tops
cup <- tr5 %>%
  dplyr::filter(all_y_f1 > .5) %>%
  dplyr::top_n(10, all_y_auc)
dcup <- tops5_to_upset(cup)

cd5all <- qs::qread(fs::path(outdir, "cd4_mod5.qs")) %>% # Long-running
  dplyr::filter(model_name %in% dcup$model_name) %>%
  dplyr::select(model_name, all_cd4_auc = all_y_auc) %>%
  unique()
dcupa <- dcup %>%
  dplyr::left_join(cd5all, by = c("model_name")) %>%
  tidyr::pivot_longer(c(all_y_auc, all_cd4_auc))

# Exploratory upset plot
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

# Recursive feature elimination (RFE) exploratory plots
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

# Exploratory: Dimension reduction on RFE results
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

# Panel B - Pseudobulk Expression
cd4data <- qread(here::here("gdm/results/cd4_model_dataset.qs"))[[2]]
c4plot <- ggplot(cd4data, aes(x = RPS18, y = `MT-CO3`)) +
  geom_point(aes(color = factor(status))) +
  theme_bw() +
  scale_color_manual(values = as.character(gdm_pal[c("orange", "black")]), labels = c("GDM", "Control")) +
  labs(color = "") +
  theme(legend.position = "none")
ggsave(fs::path(figdir, "c4_top2genes.png"), c4plot)

# Panel C - PCA on Wang and STIRM
mkl <- qs::qread(fs::path(outdir, "markers.qs"))$markers
wangmark <- load_wang_bulk("gdm/data/GSE154414_Expression_Gene.csv") %>% filter(gene %in% mkl)
wamat <- tidyr::pivot_wider()
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

# CD4 AUC scores (Described in text only)
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
# Exploratory figure
p_cmod <- ggplot(cmod_ranks, aes(y = model_name, x = mauc)) +
  geom_col() +
  theme_bw()
ggsave(fs::path(figdir, "cd4_smod_auc.png"), p_cmod)

# Single marker AUC and F1 scores
stirms_smod <- qs::qread(fs::path(outdir, "stirms_smod.qs"))
wang_smod <- qs::qread(fs::path(outdir, "wang_smod.qs"))
afscores <- rbind(
  dplyr::select(stirms_smod, model_name, data_name, all_y_auc, all_y_f1),
  dplyr::select(wang_smod, model_name, data_name, all_y_auc, all_y_f1)
) %>%
  unique() %>%
  dplyr::mutate(model_name = str_remove(model_name, "status;"))
afscores$datan2 <- ifelse(afscores$data_name == "stirms_wide", "Stirm et al.", "Wang et al.")

# Cache for supplementary table
data.table::fwrite(dplyr::filter(afscores, all_y_auc >= 0.7), "gdm/results/stirm_wang_top_supplement.csv")
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

# Raw expression values
txtdf <- data.table::fread(fs::path("gdm/results/stirm_wang_top_supplement.csv"))

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

## Figures ----

# Cache figure plot data for downstream use
figcache <- list(
  p_odds,
  p_cd4_upset,
  prd_roc_plot,
  prd_f1_plots,
  p_rfe,
  p_dimreg,
  c4plot,
  afplot,
  p_expr2
)
# Note: p_expr2 not included in figcache
qs::qsave(figcache, fs::path(outdir, "gdm_covar_figcache_v3.qs"))

# Create Figure panel

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



## Export ----

figcache <- qs::qread(fs::path(outdir, "gdm_covar_figcache_v2.qs"))

p_odds <- figcache[[1]]
data.table::fwrite(p_odds$data, fs::path(zdir, "f4_odds.csv"))

c4plot <- figcache[[7]]
data.table::fwrite(
  c4plot$data %>% dplyr::select(id, status, RPS18, `MT-CO3`),
  fs::path(zdir, "f4_genes_pair.csv")
)

afplot <- figcache[[8]]
names(afplot$data) <- c("gene_model", "data", "auc", "f1", "data_name")
data.table::fwrite(afplot$data, fs::path(zdir, "f4_auc_f1.csv"))

expr <- figcache[[9]]
data.table::fwrite(expr$data, fs::path(zdir, "f4_expression.csv"))
