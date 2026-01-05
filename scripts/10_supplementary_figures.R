# Create supplementary figures and export underlying data objects

# Propeller testing ----

# Library ----
source("code/plots.R")
source("code/subtyping_de.R")
library(data.table)

# Variables ----
outdir <- "gdm/results"
datadir <- "gdm/data"
figdir <- fs::dir_create("gdm/figures/main")
zdir <- fs::dir_create("gdm/export/") # Exports datasets underlying figures

abfigobs <- qs::qread("gdm/results/ab_test_fig.qs")
abcd4 <- abfigobs$cd4_ab

f1pal_order <- c("darkpurp", "orange", "darkred", "blue", "ogreen", "darkgrey", "purple", "pink", "darkgrey", "teal")
f2pal <- gdm_dimreg_pal[match(f1pal_order, names_gdm_dimreg_pal)]

abcd4_2 <- abcd4[[2]] + scale_fill_manual(values = f2pal) +
  labs(fill = "scRNA\ncluster", title = "CD4+ T\nscRNA") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45))

cd4patplot <- abcd4[[4]] + labs(x = "Patient", fill = "Disease\nStatus")
cd4abplt <- abcd4[[3]] + labs(x = "Ratio\n(GDM:CONTROL)")

cd4_f2 <- patchwork::wrap_plots(
  abcd4_2,
  # patchwork::free(cd4patplot, type='label'),
  patchwork::free(cd4abplt, type = "label"),
  nrow = 1,
  widths = c(0.2, 0.4)
) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") & theme(
  legend.position = "bottom",
  legend.direction = "horizontal", # Make legend horizontal
  legend.title = element_text(hjust = 0.5), # Center title (optional)
  legend.box = "vertical", # Stack title above legend
  plot.margin = margin(5, 5, 5, 5)
)
ggsave(fs::path(figdir, "Figure_S_CD4abundance.png"), cd4_f2, width = 8.4, height = 6)


# CD4+ T cells ----

# Compose Supplementary Figure 1: CD4+ T cells
source("code/plots.R")
source("code/subtyping_de.R")
library(qs)

outdir <- here::here("gdm/results")
figdir <- fs::dir_create("gdm/figures/main")
gdm <- qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 <- gdm$cd4
treg <- gdm$treg

# Plot TSNE for CD4+ T
Idents(cd4) <- "cluster"
cd4_tsne <- show_tsneplot(cd4, pal = cd4_dimreg_pal) + labs(title = "CD4+ T Cells")

# Show highlight features
cd4_highlight <- c("CCR7", "IL2RA", "CCR6", "CXCR3", "S100A4", "FOXP3")
p_cd4high <- FeaturePlot(cd4, features = cd4_highlight, cols = c("lightgray", gdm_pal[["orange"]]), reduction = "tsne") + NoLegend()

# Show heatmap for CD4
top_markers <- readr::read_csv(fs::path("gdm/results/", "top_subtyping_de_markers.csv"))
markers <- de_top_n_genes(dplyr::filter(top_markers, label == "cd4"), 2)
markers <- c(markers, "FOXP3", cd4_highlight)
cd4_avc <- qs::qread(fs::path(outdir, "gdm_average_expr.qs"))$cd4_avc
cd4_hm <- subtype_heatmap(cd4_avc, markers, is.treg = F)

# Humanin+ boxplots
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

pp <- qs::qread(fs::path(outdir, "propeller.qs"))
tpropeller <- pp[[1]]
cpropeller <- pp[[2]]
p_ab_test_cd4 <- plot_propeller_res(cpropeller %>% filter(cluster == "Memory-1"), "CD4", gdm_pal[["pink"]])
p_ab_test_cd4

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


# GDM analysis marker genes: Supplementary Table ----

# Get a unique list of marker genes
mkl <- qs::qread(fs::path("gdm/results", "markers.qs"))
markers <- mkl$markers
marklist <- mkl$markerlist

# For marker gene sources, setup a look up table for their supplement name
codes <- c("rdeg_treg", "rdeg_cd4", "nact_treg", "hub_string", "humanin")
stopifnot(all(codes %in% names(marklist)))
code_rename <- c("Treg DEG", "CD4+ DEG", "Treg Naive-Act-2 DEG", "Treg Naive-Act-2 Hub Gene", "Treg Humanin+ DEG")
names(code_rename) <- codes

# Get marker gene sources per gene
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

# Create marker gene results table
mgt <- data.frame(
  "Gene" = markers,
  "Criteria" = sapply(markers, function(e) gmgspg(e))
)
data.table::fwrite(mgt, "gdm/results_v2/marker_genes.csv")

# Humanin scores ----
library(Seurat)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(dplyr)

seu <- qs::qread(fs::path("gdm/results/", "gdm_subtyped.qs"))
treg <- seu$treg

treg$hvsall <- ifelse(treg$cluster == "Humanin+", "Humanin+", "Other")
treg$S.Score_1 <- treg$S.Score
treg$G2M.Score_1 <- treg$G2M.Score

# Apply again
treg <- CellCycleScoring(treg,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes
)

metadata <- treg@meta.data

# Create new variable with GDM features
metadata$xval <- paste0(metadata$cond_cln, "\n", metadata$hvsall)

# Calculate p-values
p1a <- wilcox.test(
  x = dplyr::filter(metadata, hvsall == "Humanin+") %>% pull(S.Score),
  y = dplyr::filter(metadata, hvsall == "Other") %>% pull(S.Score),
  paired = FALSE
)
p1b <- wilcox.test(
  x = dplyr::filter(metadata, hvsall == "Humanin+") %>% pull(S.Score_1),
  y = dplyr::filter(metadata, hvsall == "Other") %>% pull(S.Score_1),
  paired = FALSE
)
p2a <- wilcox.test(
  x = dplyr::filter(metadata, hvsall == "Humanin+") %>% pull(G2M.Score),
  y = dplyr::filter(metadata, hvsall == "Other") %>% pull(G2M.Score),
  paired = FALSE
)
p2b <- wilcox.test(
  x = dplyr::filter(metadata, hvsall == "Humanin+") %>% pull(G2M.Score_1),
  y = dplyr::filter(metadata, hvsall == "Other") %>% pull(G2M.Score_1),
  paired = FALSE
)
padj <- p.adjust(c(p1b$p.value, p1a$p.value, p2b$p.value, p2a$p.value), method = "fdr")
padj <- format(padj, digits = 2)

plot_ccscores <- function(field = "S.Score", pannotations = padj[1], title = "CellCycle:S") {
  ggplot(metadata, aes_string(x = "hvsall", y = field, fill = "hvsall")) +
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    theme_classic() +
    scale_fill_manual(values = c(scales::alpha("#004949", 0.5), "lightgrey")) +
    labs(
      x = "Cluster", y = field,
      title = title, fill = ""
    ) +
    theme(legend.position = "none") +
    geom_signif(
      comparisons = list(
        c("Humanin+", "Other")
      ), y_position = c(0.3),
      annotations = pannotations
    ) +
    scale_y_continuous(limits = c(0, 0.35))
}

p_sbefore <- plot_ccscores("S.Score_1", padj[1], "S.Score Before")
p_safter <- plot_ccscores("S.Score", padj[2], "S.Score After")
p_mbefore <- plot_ccscores("G2M.Score_1", padj[3], "G2M.Score Before")
p_mafter <- plot_ccscores("G2M.Score", padj[4], "G2M.Score After")

pfig <- (p_sbefore + p_safter) / (p_mbefore + p_mafter)
ggsave("data/hu_alignment/humanin_scores.png", pfig, width = 8, height = 6, dpi = 300)

# Run VDJ clonotype figures ----

# Library ----
library(Seurat)
library(scRepertoire)
library(dplyr)
library(ggplot2)
library(data.table)
library(parallel)
source("code/preprocessing_qc.R")
source("code/vdj_helpers.R")

# Variables ----
indir <- fs::dir_create(here::here("data/seurat/220802_Shangaris"))
outdir <- fs::dir_create(here::here("gdm/results/vdj/"))
figdir <- fs::dir_create(here::here("gdm/figures/vdj/"))
meta_csv <- here::here("data/metadata/220721_metadata_cln.csv")
tenx_file <- fs::path(indir, "seurat_tenx_clustered.qs")
tcr_file <- fs::path(indir, "screp_combined_tcr.qs")

# Main ----

# Load scRNA data
tenx <- qs::qread(tenx_file)
DefaultAssay(tenx) <- "RNA"
tenx <- AddRutepoMeta(tenx, meta_csv, cols_regex = "developed|condition|control")
cl <- qs::qread(fs::path(indir, "screp_combined_tcr.qs")) # clonotypes - source: ./legacy/scripts/220907_VDJ_clonotypes.R
gdm <- qs::qread(fs::path("gdm/results/", "gdm_subtyped.qs"))
cd4 <- gdm$cd4
treg <- gdm$treg

# Select VDJ data from samples used for scRNA analysis
keep_samples <- unique(unique(c(treg$multi_q, cd4$multi_q)))
has_keep <- function(x) {
  checked_names <- sapply(keep_samples, function(i) {
    x_clean <- stringr::str_remove(x, "-Treg")
    x_clean <- stringr::str_remove(x_clean, "-CD4")
    stringr::str_detect(x_clean, i)
  })
  any(checked_names)
}
keep_vector <- sapply(names(cl), has_keep)
cl <- cl[keep_vector]

# Order TCR data
ord <- data.table(samples = names(cl))
mng <- tidyr::separate(ord, samples, c("patient", "condition"), sep = "_", remove = F, extra = "merge") %>%
  tidyr::separate(patient, into = c("run", "hashtag"), sep = -5, remove = F) %>%
  dplyr::mutate(run = str_remove(run, "-$")) %>%
  tidyr::separate(run, into = c("base", "ix", "ctype"), sep = "-", remove = F) %>%
  dplyr::arrange(ctype, condition, base, hashtag) %>%
  dplyr::mutate(cond_cln = ifelse(condition == "CONTROL", "CONTROL", "GDM"))
cl <- cl[mng$samples]

# Relative abundance analysis per patient
p_homeo <- patchwork::wrap_plots(lapply(conds, plt_homeo, cl = cl), ncol = 1)
ggsave(fs::path(figdir, "vdj_homestasis.png"), p_homeo, width = 8, height = 8)

# Relative abaundance per patient per group
data <- do.call(rbind, lapply(p_homeo, function(x) x$data))
names(data) <- c("samples", "Clonotype Group", "Relative Abundance")
data <- data %>% left_join(mng)
p_cd4 <- ggplot(
  data %>% dplyr::filter(ctype == "CD4"),
  aes(x = patient, y = `Relative Abundance`, fill = `Clonotype Group`)
) +
  geom_col() +
  facet_wrap(~cond_cln, scales = "free") +
  labs(title = "CD4", x = "Patient") +
  scale_fill_viridis_d(option = "plasma") +
  coord_flip()
p_treg <- ggplot(
  data %>% dplyr::filter(ctype == "Treg"),
  aes(x = patient, y = `Relative Abundance`, fill = `Clonotype Group`)
) +
  geom_col() +
  facet_wrap(~cond_cln, scales = "free") +
  labs(title = "Treg", x = "Patient") +
  scale_fill_viridis_d(option = "plasma") +
  coord_flip()
ggsave(fs::path(figdir, "vdj_homestasis_clean.png"), p_treg / p_cd4, width = 12, height = 8)

# Relative abundance per group (GDM vs Treg)
clcomb <- do.call(rbind, cl)
mng$sample <- mng$patient
clcomb <- clcomb %>% left_join(dplyr::rename(mng, "values" = "samples"))
f_tregs <- clonalHomeostasis(clcomb %>% dplyr::filter(ctype == "Treg"),
  cloneCall = "gene",
  cloneTypes = c(
    Rare = 1e-04,
    Small = 0.001,
    Medium = 0.01,
    Large = 0.1,
    Hyperexpanded = 1
  ),
  group.by = "cond_cln"
) +
  labs(title = "Treg") + coord_flip()
f_cd4 <- clonalHomeostasis(clcomb %>% dplyr::filter(ctype == "CD4"),
  cloneCall = "gene",
  cloneTypes = c(
    Rare = 1e-04,
    Small = 0.001,
    Medium = 0.01,
    Large = 0.1,
    Hyperexpanded = 1
  ),
  group.by = "cond_cln"
) +
  labs(title = "CD4") + coord_flip()
p_homeo <- patchwork::wrap_plots(list(f_tregs, f_cd4), ncol = 1)
ggsave(fs::path(figdir, "vdj_homestasis_group.png"), p_homeo, width = 8, height = 8)

# Visualise the Humanin Serum Elisa tests ----

# Lirbary
library(data.table)
library(ggplot2)
library(janitor)
library(dplyr)
library(patchwork)
library(ggsignif)
source("code/plots.R")

# Load data
el <- fread("gdm/data/HUMANIN_ELIZA_CLEAN.csv")
el <- janitor::clean_names(el)
el <- dplyr::filter(el, source != "")
el <- dplyr::mutate(el,
  cond_cln = ifelse(condition == "Healthy", "Control", condition),
  period = ifelse(source == "growth", "36 Weeks", "12 Weeks")
)

# Visualise
p0 <- ggplot(el %>% dplyr::filter(period == "12 Weeks"), aes(y = humanin_ng_m_l, x = cond_cln, color = cond_cln)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.8, width = 0.15) +
  theme_bw() +
  scale_color_manual(values = c("grey", gdm_pal[["orange"]])) +
  labs(y = "Serum Humanin (ng/ml)", title = "12 Weeks", x = "") +
  theme(legend.position = "none") +
  geom_signif(
    comparisons = list(c("GDM", "Control")),
    test = "t.test",
    test.args = list(alternative = "greater"),
    color = "black"
  )
p1 <- ggplot(el %>% dplyr::filter(period == "36 Weeks"), aes(y = humanin_ng_m_l, x = cond_cln, color = cond_cln)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.8, width = 0.15) +
  theme_bw() +
  scale_color_manual(values = c("grey", gdm_pal[["orange"]])) +
  labs(y = "Serum Humanin (ng/ml)", title = "36 Weeks", x = "") +
  theme(legend.position = "none") +
  geom_signif(
    comparisons = list(c("GDM", "Control")),
    test = "t.test",
    test.args = list(alternative = "greater"),
    color = "black"
  )
plt <- p0 + p1

ggsave("humanin_elisa_vis.png", plt, width = 8, height = 5, dpi = 300)

# TXNIP expression in single cell datasets ----
gdm <- qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 <- gdm$cd4
treg <- gdm$treg
rm(gdm)
Idents(treg) <- "cluster"
treg_txnip <- VlnPlot(treg, "TXNIP", split.by = "control", log = F, add.noise = T, slot = "scale.data", pt.size = 0.00005, cols = c("orange", "lightgrey")) + labs(title = "Treg | TXNIP")
Idents(cd4) <- "cluster"
cd4_txnip <- VlnPlot(cd4, "TXNIP", split.by = "control", log = F, add.noise = T, slot = "scale.data", pt.size = 0.00005, cols = c("orange", "lightgrey")) + labs(title = "CD4+ | TXNIP")
cpatch <- patchwork::wrap_plots(treg_txnip, cd4_txnip, ncol = 1) +
  plot_annotation(tag_levels = "A")
ggsave(filename = fs::path(figdir, "TXNIP.png"), plot = cpatch, width = 8, height = 8)
