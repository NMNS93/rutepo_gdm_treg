source("code/glm_meta_model.R")
source("code/plots.R")
source("code/stringdb.R")
library(logger)
library(MASS)
library(dplyr)
library(logger)

# ==============================================================================

# Setup ----
outdir <- here::here("gdm/results_v3")
figdir <- fs::dir_create("gdm/figures_v3/cellpheondb")
gdm <- qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 <- gdm$cd4
treg <- gdm$treg
rm(gdm)
seurat_obj <- subset(cd4, subset = (cluster != "MALAT1+"))

# Split the Seurat object by 'control' field (CASE and CONTROL)
seurat_case <- subset(seurat_obj, subset = control == "CASE")
seurat_control <- subset(seurat_obj, subset = control == "CONTROL")

# ---- For CASE ----
# Get the normalized counts for CASE
expr_matrix_case <- as.data.frame(GetAssayData(seurat_case, slot = "data"))
write.table(expr_matrix_case, file = "data/cellphonedb/cd4_case_expression_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Extract cell type information for CASE
metadata_case <- data.frame(
    Cell = colnames(seurat_case),
    cell_type = seurat_case$cluster
) # Replace 'cluster' with the actual column name for cell type

# Save the metadata file for CASE
write.table(metadata_case, file = "data/cellphonedb/cd4_case_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# ---- For CONTROL ----
# Get the normalized counts for CONTROL
expr_matrix_control <- as.data.frame(GetAssayData(seurat_control, slot = "data"))
write.table(expr_matrix_control, file = "data/cellphonedb/cd4_control_expression_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Extract cell type information for CONTROL
metadata_control <- data.frame(
    Cell = colnames(seurat_control),
    cell_type = seurat_control$cluster
) # Replace 'cluster' with the actual column name for cell type

# Save the metadata file for CONTROL
write.table(metadata_control, file = "data/cellphonedb/cd4_control_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Setup DE file - This will be the same in both analyses. It identifies genes DE in the GDM FOXP3+ Tregs
## Define DEGs as those significantly altered in a Treg subset
tde <- fread("gdm/results_v3/treg_within_cluster_de.csv.gz")
tdata <- tde[stringr::str_detect(cluster, "Effector") & is_signif == T]
tdata$cluster <- "FOXP3+"
tdata <- tdata[, .(cluster, gene)] %>% unique()
data.table::fwrite(tdata, "data/cellphonedb/cd4_degs.txt", sep = "\t")

# Pause --- Here we manually run cellphonedb.py separately.

# CellphoneDB plots ----
library(ktplots)
pvals <- read.delim("gdm/results_v3/cellphonedb/case/degs_analysis_relevant_interactions_09_28_2024_140844.txt")
means <- read.delim("gdm/results_v3/cellphonedb/case/degs_analysis_means_09_28_2024_140844.txt")
decon <- read.delim("gdm/results_v3/cellphonedb/case/degs_analysis_deconvoluted_09_28_2024_140844.txt")
obj <- seurat_case

plot_cpdb_heatmap(pvals = pvals, cellheight = 10, cellwidth = 10)
plot_cpdb_heatmap(pvals = pvals, cellheight = 10, cellwidth = 10, deg_analysis = T)

# Control
pvals <- read.delim("gdm/results_v3/cellphonedb/control/degs_analysis_relevant_interactions_09_28_2024_142839.txt")
pvals <- read.delim('gdm/results_v3/cellphonedb/control/degs_analysis_interaction_scores_09_28_2024_142839.txt')
means <- read.delim("gdm/results_v3/cellphonedb/control/degs_analysis_means_09_28_2024_142839.txt")
decon <- read.delim("gdm/results_v3/cellphonedb/control/degs_analysis_deconvoluted_09_28_2024_142839.txt")
obj <- seurat_control

plot_cpdb_heatmap(pvals = pvals, cellheight = 10, cellwidth = 10, deg_analysis = F)
plot_cpdb(
    scdata = seurat_control,
    cell_type1 = "FOXP3",
    cell_type2 = ".", # this means all cell-types
    celltype_key = "cluster",
    means = means,
    pvals = pvals,
    # genes=c("PTPRC", "TNFSF13"),
    title = "Interactions",
)

# Stat CASE
pvals <- read.delim('gdm/results_v3/cellphonedb_stat/statistical_analysis_pvalues_09_28_2024_182652.txt')
rel <- fread('gdm/results_v3/cellphonedb/case/degs_analysis_relevant_interactions_09_30_2024_025638.txt')
pints <- fread("gdm/results_v3/cellphonedb/case/degs_analysis_interaction_scores_09_30_2024_025638.txt")
sints <- fread("gdm/results_v3/cellphonedb/case/degs_analysis_deconvoluted_percents_09_30_2024_025638.txt")
smeans <- fread("gdm/results_v3/cellphonedb/case/degs_analysis_significant_means_09_30_2024_025638.txt")
means <- fread("gdm/results_v3/cellphonedb/case/degs_analysis_means_09_30_2024_025638.txt")
decon <- fread("gdm/results_v3/cellphonedb/case/degs_analysis_deconvoluted_09_30_2024_025638.txt")
obj <- seurat_case


# As no hits were found, save the cpdb results for FOXP3+ cells.
res <- dplyr::select(means, c(1:13), contains("FOXP3")) %>%
    mutate(dataset = "CD4+_GDM") %>%
    tidyr::pivot_longer(cols = contains("FOXP3"), names_to = "pairs", values_to = "means") %>%
    data.table()
issig <- dplyr::select(smeans, c(1:13), contains("FOXP3")) %>%
    mutate(dataset = "CD4+_GDM") %>%
    tidyr::pivot_longer(cols = contains("FOXP3"), names_to = "pairs", values_to = "issig") %>%
    dplyr::filter(!is.na(issig))
res <- merge(res, issig, all.x = T)
rlints <- dplyr::select(rel, c(1:13), contains("FOXP3")) %>%
    mutate(dataset = "CD4+_GDM") %>%
    tidyr::pivot_longer(cols = contains("FOXP3"), names_to = "pairs", values_to = "relevant") %>%
    data.table() %>%
    dplyr::filter(relevant > 0) %>%
    pull(id_cp_interaction)

# Show relevant (if any)
gp <- ggplot(drel, aes(x = pairs, y = interacting_pair)) +
    geom_point(aes(size = means, fill = sig_bool), shape = 21, color = "black") +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
ggsave("test.png", gp)

# Show interactions for FOXP3+
res <- dplyr::select(means, c(1:13), contains("FOXP3")) %>%
    tidyr::pivot_longer(cols = contains("FOXP3"), names_to = "pairs", values_to = "means") %>%
    data.table()
gp <- ggplot(res, aes(x = pairs, y = interacting_pair)) +
    geom_point(aes(size = means, fill = means), shape = 21, color = "black") +
    # scale_fill_brewer(palette='Set2') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
ggsave("ints.png", gp, width = 10, height = 10)
data.table::fwrite(res, "gdm/results_v3/case_cellphonedb.csv")
