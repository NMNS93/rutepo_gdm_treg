# Define a predictive model for GDM using RNA-seq data
# Combines metadata, pseduobulk scRNA and bulk blood mRNA

source("code/glm_meta_model.R")
source("code/plots.R")
source("code/stringdb.R")
library(logger)
library(MASS)
library(dplyr)
library(logger)

# ==============================================================================

# Setup ----
outdir = here::here("gdm/results_v3")
figdir = fs::dir_create("gdm/figures_v3/covariates")
gdm = qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 = gdm$cd4; treg = gdm$treg; rm(gdm)

# Metadata ----
log_info("Setting up Metadata")

# Extract GDM risk factors from metadata to use as predictors

# Load and format metadata
meta_csv = here::here("data/metadata/220721_metadata_cln.csv")
meta = readr::read_csv(meta_csv)
meta <- blood_pressure_refactor(meta)
meta <- medication_reformat(meta)
meta <- dplyr::select(meta, meta_select_vars)

# Drop sample without BMI or MAP data
exclude_sample = meta$multi_q[!(complete.cases(meta))]
log_info("Excluding sample {exclude_sample}")
meta = dplyr::filter(meta, multi_q!=exclude_sample)

# Add case-control status
meta <- add_case_control_seurat(meta, cd4, treg)

# Fix metadata names
meta <- dplyr::rename(meta, id=multi_q, status=cond_cln)
meta$status_c = meta$status
meta$status = ifelse(meta$status=='CONTROL', 0, 1)
  
# Save and plot data
p_meta <- plot_gdm_meta(meta)
ggsave(fs::path(figdir, "meta_data_summary.png"), p_meta$plot, width=5, height=10)
data.table::fwrite(meta, fs::path(outdir, "meta_data.csv"))
meta <- fread(fs::path(outdir, "meta_data.csv"))

# Treg subset marker genes ----
log_info("Treg subset markers")
# Define DEGs as either Treg subtype marker genes or Treg within-cluster DEGs
markerlist <- list()

# Load differential expressed genes
top_markers = readr::read_csv(fs::path(outdir, "top_subtyping_de_markers.csv")) %>% dplyr::filter(label=='treg')
de = qs::qread(fs::path(outdir, "we_treg_cln.qs"))

# List markers
markerlist = list()
markerlist$clustmarks <- top_markers %>% 
  filter(
    CONTROL_p_val_adj <= 0.01 & CASE_p_val_adj <= 0.01 & 
    abs(CONTROL_avg_log2FC) >= 0.5 & abs(CASE_avg_log2FC)>=.5
    ) %>% 
  pull(gene) %>% unique()

markerlist$de_marks <- de$wta %>% filter(is_signif==T) %>% pull(gene) %>% unique()

# List unique markers
markers_rna = unique(unlist(markerlist))
markers = gsub('AL138963\\.4', 'TPT1', markers_rna)

# Cache
qs::qsave(
  list('markerlist'=markerlist, 'markers'=markers),
  fs::path(outdir, "markers.qs")
)
mkl = qs::qread(fs::path(outdir, "markers.qs"))
markers=mkl$markers

# Create CD4+ and Treg pseudobulks ----
log_info("Creating CD4 pseudobulks")
avcd4dt <- make_scrna_model_dataset(cd4, markers_rna)
names(avcd4dt) = gsub('AL138963\\.4', 'TPT1', names(avcd4dt))
cd4_dataset <- merge(avcd4dt, meta[, c("status", "id")], all.x=T)
names(cd4_dataset) = gsub('\\.', '-', names(cd4_dataset))
qs::qsave( list(avcd4dt, cd4_dataset), fs::path(outdir, "cd4_model_dataset.qs") )
ll = qread(fs::path(outdir, "cd4_model_dataset.qs"))
avcd4dt = ll[[1]]; cd4_dataset = ll[[2]]

avtregdt <- make_scrna_model_dataset(treg, markers_rna)
names(avtregdt) = gsub('AL138963\\.4', 'TPT1', names(avtregdt))
treg_dataset <- merge(avtregdt, meta[, c("status", "id")], all.x=T)
names(treg_dataset) = gsub('\\.', '-', names(treg_dataset))
qs::qsave( list(avtregdt, treg_dataset), fs::path(outdir, "treg_model_dataset.qs") )
mm = qread(fs::path(outdir, "treg_model_dataset.qs"))
avtregdt = mm[[1]]; treg_dataset = mm[[2]]

# Load Stirm et al. (2018) bulk RNA-seq data ----
log_info("Loading bulk stirm")
# Dataset of bulk RNA-seq from patients with and without GDM.
# Source: https://www.nature.com/articles/s41598-018-19200-9#Sec11
stirm <- load_stirm_bulk("data/mrna/GSE92772_expression_mrna.tsv", "data/mrna/GSE92772_series_matrix.txt")
stirms = dplyr::filter(stirm, gene %in% markers)
stirms_wide = tidyr::pivot_wider(stirms, names_from="gene", values_from="expr") %>%
    dplyr::mutate(status = ifelse(status=="GDM", 1, 0))
qs::qsave(stirms_wide, fs::path(outdir, 'strims_wide.qs')) # For awards

# Generate exploratory plots showing marker genes in stirms dataset
pbox = plot_stirms_box(stirms, figdir)
pgene = p_wrap_stirms_gene(stirms, figdir)

# Load Wang et al. (2021) bulk RNA-seq data ----
log_info("Loading bulk Wang")

# Source: https://www.sciencedirect.com/science/article/pii/S088875432100077X
wangl = load_wang_bulk("gdm/data/GSE154414_Expression_Gene.csv")
wangls = dplyr::filter(wangl, gene %in% markers)
wang_wide = tidyr::pivot_wider(wangls, names_from='gene', values_from='expression')
qs::qsave(wang_wide, fs::path(outdir, 'wang_wide.qs'))

# Modeling Background ----

# We aim to perform LOO-CV, measuring AUC scores, AIC, classification accuracy on:
# 1) The baseline model with metadata only
# 2) The average CD4 dataset with combinations of marker genes
# 3) The average Treg dataset with combinations of marker genes
# 4) The top 10 best-performing models in the average CD4 dataset (#2)

# We define a pipeline that accepts a dataset, a model object defining which fields to fit.
# The pipeline will always predict on the 'status' field. It returns an object carrying LOO-CV scores.
# The object is accessible with the attribs: model, auc, aic, classification
# Note that we test combinations of 1, 3 and 5 genes

# Model metadata ----
log_info("Fitting model metadata")

# Select variables. Treatment such as metformin are excluded as they are given after GDM determined.
independent_vars <- c("ethnicity", "bmi", "map")
# Generate all combinations of variables with target 'status' included
var_selectors <- build_variable_selectors(independent_vars)

# Fit models and collect metrics
meta_mod = fit_and_collect_model_metrics(meta, "meta", var_selectors, merge=T)
qs::qsave(meta_mod, fs::path(outdir, "meta_mod.qs"))
meta_mod <- qs::qread(fs::path(outdir, 'meta_mod.qs'))

# Calculate metadata odds ratios ----
log_info("Calculating odds ratios")
odds_data <- metadata_odds_ratio(meta)
data.table::fwrite(odds_data, fs::path(outdir, 'meta_odds.csv'))

# Single-gene models (bulk and pseudobulk) ----

# Stirm et al
log_info("Running single-gene stirm")
var_sel_stirms_single = build_variable_selectors(unique(stirms$gene), sizes=1)
stirms_smod = fit_and_collect_model_metrics(stirms_wide, "stirms_wide", var_sel_stirms_single, merge=T)
qs::qsave(stirms_smod, fs::path(outdir, "stirms_smod.qs"))
stirms_smod <- qs::qread(fs::path(outdir, "stirms_smod.qs"))

# Wang et al
log_info("Running single-gene wang")
var_sel_wang = build_variable_selectors(unique(wangls$gene), sizes=1)
wang_smod = fit_and_collect_model_metrics(wang_wide, "wang_wide", var_sel_wang, merge=T)
qs::qsave(wang_smod, fs::path(outdir, "wang_smod.qs"))
wang_smod = qs::qread(fs::path(outdir, "wang_smod.qs"))

# CD4 Pseu
log_info("Running single-gene cd4")
var_sel_cd4 = build_variable_selectors(setdiff(names(cd4_dataset), c("id", "status")), sizes=1)
cd4_smod =  fit_and_collect_model_metrics(cd4_dataset, "cd4_wide", var_sel_cd4, merge=T)
qs::qsave(cd4_smod, fs::path(outdir, "cd4_smod.qs"))
cd4_smod <- qs::qread(fs::path(outdir, "cd4_smod.qs"))

# Treg Pseu
log_info("Running single-gene treg")
var_sel_treg = build_variable_selectors(setdiff(names(treg_dataset), c("id", "status")), sizes=1)
treg_smod =  fit_and_collect_model_metrics(treg_dataset, "treg_wide", var_sel_treg, merge=T)
qs::qsave(treg_smod, fs::path(outdir, "treg_smod.qs"))
treg_smod <- qs::qread(fs::path(outdir, "treg_smod.qs"))

# Re-run using the top across all samples (for plotting downstream)
roc_data_ss <- smod_to_roc_ss(wang_smod, stirms_smod, cd4_smod, treg_smod)
qs::qsave(roc_data_ss, fs::path(outdir, "roc_data_ss_cache.qs"))
nmodname = roc_data_ss$model_name %>% unique %>% gsub('\\.', '-', .)
var_sel_sstops <- build_variable_selectors(nmodname, sizes=1)

smod_all_tops <- list(
  fit_and_collect_model_metrics(stirms_wide, 'stirms_wide', var_sel_sstops, merge=T),
  fit_and_collect_model_metrics(wang_wide, 'wang_wide', var_sel_sstops, merge=T),
  fit_and_collect_model_metrics(cd4_dataset, 'cd4_wide', var_sel_sstops, merge=T),
  fit_and_collect_model_metrics(treg_dataset, 'treg_wide', var_sel_sstops, merge=T)
)
qs::qsave(smod_all_tops, fs::path(outdir, "smod_all_tops.qs"))

roc_data_ss2 <- smod_to_roc_ss(
  smod_all_tops[[2]], smod_all_tops[[1]], smod_all_tops[[3]], smod_all_tops[[4]]
)
qs::qsave(roc_data_ss2, fs::path(outdir, "roc_data_ss2_cache.qs"))

## Recursive Feature Elimination ----
log_info("RFE pipes")

rfe_pipe <- function(x){
  X = x %>% dplyr::select(-id, -status) %>% as.matrix()
  Y = x %>% dplyr::select(status) %>% as.matrix()
  rfctrl <- rfeControl(functions = rfFuncs, method = "LOOCV")
  
  # Run RFE with Random Forest (Note Logistic Regression fails)
  rferf <- rfe(X, Y, sizes = 1:5, rfeControl = rfctrl)
  
  # Get AUC
  rfauc <- ModelMetrics::auc(Y[,1], predict(rferf$fit, X))
  log_info('RF AUC {rfauc} | top: {paste0(rferf$optVariables, collapse=";")}')
  return(
    list(rferf, rfauc)
  )
}

log_info("RFE stirm")
stirm_rfe = rfe_pipe(stirms_wide)
log_info("RFE wang")
wang_rfe = rfe_pipe(wang_wide %>% dplyr::select(-status_c))
log_info("RFE cd4")
cd4_rfe = rfe_pipe(cd4_dataset)
log_info("RFE treg")
treg_rfe = rfe_pipe(treg_dataset)

qs::qsave(
  list('w'=wang_rfe, 's'=stirm_rfe, 'c'=cd4_rfe, 't'=treg_rfe),
  fs::path(outdir, "rfe_singles.qs")
)


# Model CD4 pseudobulk ----

# Do the top 20 genes in Treg single have more predictive power
# in CD4 pseudobulks and/or in bulk RNA-seq data?
treg_tops <- get_mod_tops(treg_smod)
treg_tops_ngenes = treg_tops$tops %>% dplyr::filter(all_y_auc >= .5) %>% 
  arrange(-all_y_auc) %>% slice(1:12) %>% pull(model_name) %>%
  stringr::str_split(., ';') %>%
  unlist() %>% purrr::discard(~stringr::str_detect(.x, 'status'))
treg_tops_varsel = build_variable_selectors(treg_tops_ngenes, sizes=5)

log_info("Modelling CD4 pseudobulk")

# Fit models and collect metrics
cd4_mod_files = fit_and_collect_model_metrics(cd4_dataset, "cd4_pseu", treg_tops_varsel)
cd4_mod5 = merge_metrics_mem_aware(cd4_mod_files, fs::path(outdir, "cd4_mod5.qs"))
cd4_mod5 = qs::qread(fs::path(outdir, "cd4_mod5.qs"))
cd4_tops5 = get_mod_tops(cd4_mod5)
qs::qsave(cd4_tops5, fs::path(outdir, "cd4_tops5.qs"))

# Model Treg pseudobulk ----
log_info("Modelling Treg pseudobulk")
# Fit models and collect metrics
treg_mod_files = fit_and_collect_model_metrics(treg_dataset, "treg_pseu", treg_tops_varsel)
treg_mod5 = merge_metrics_mem_aware(treg_mod_files, fs::path(outdir, "treg_mod5.qs"))
treg_mod5 = qs::qread(fs::path(outdir, "treg_mod5.qs"))
treg_tops5 = get_mod_tops(treg_mod5)
qs::qsave(treg_tops5, fs::path(outdir, "treg_tops5.qs"))
