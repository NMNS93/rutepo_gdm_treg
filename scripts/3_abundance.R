# Run differential abundance testing with speckle/propeller method

source("code/plots.R")
library(patchwork)
library(speckle)
library(dplyr)
library(logger)
outdir=fs::dir_create("gdm/results_v3/")
figdir=fs::dir_create("gdm/figures_v3/abundance")
gdm = qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
treg = gdm$treg; cd4=gdm$cd4; rm(gdm)
DefaultAssay(treg) <- "SCT"
DefaultAssay(cd4) <- "SCT"

# Propeller
# https://academic.oup.com/bioinformatics/article/38/20/4720/6675456
# https://github.com/phipsonlab/speckle
log_info("Running Treg propeller")
tpropeller <- propeller(clusters=treg$cluster, sample=treg$multi_q, group=treg$cond_cln) %>%
    dplyr::mutate(cluster=rownames(.)) %>% dplyr::mutate(PropRatio=1+(1-PropRatio)) # Flip for GDM
log_info("Completed Treg propeller")
log_info("Running CD4 propeller")
cpropeller <- propeller(clusters=cd4$cluster, sample=cd4$multi_q, group=cd4$cond_cln) %>%
  dplyr::mutate(cluster=rownames(.)) %>% dplyr::mutate(PropRatio=1+(1-PropRatio))
log_info("Completed Treg propeller")

log_info("Saving propeller results")
qs::qsave(list(tpropeller, cpropeller), fs::path(outdir, "propeller.qs"))

# Load properller results
pp = qs::qread(fs::path(outdir, 'propeller.qs'))
tpropeller = pp[[1]]
cpropeller = pp[[2]]
rm(pp)

# Plot 1: Fill by condition
log_info("Generating plots")
p_ab_fbc_treg <- plot_abundance_fbc(treg, "", pal=treg_dimreg_pal)
p_ab_fbc_cd4 <- plot_abundance_fbc(cd4, "", pal=cd4_dimreg_pal)

# Plot 2: Propeller fold-change results
p_ab_test_treg <- plot_propeller_res(tpropeller, "Treg", gdm_pal[["orange"]])
p_ab_test_cd4 <- plot_propeller_res(cpropeller, "CD4", gdm_pal[["pink"]])

# Plot 3: Patient cluster proportions
p_ab_prop_treg <- plot_pat_proportion(treg, gdm_pal[["orange"]]) + theme(legend.position = "bottom")
p_ab_prop_cd4 <- plot_pat_proportion(cd4, gdm_pal[["pink"]]) + theme(legend.position="bottom")

# Plot 4: Clusters by disease status
p_cond_treg <- plot_cond_tsne(treg, gdm_dimreg_pal[8]) + theme(legend.position="bottom")
p_cond_cd4 <- plot_cond_tsne(cd4, gdm_dimreg_pal[5])+ theme(legend.position="bottom")

# Save images to object
log_info("Saving figure datasets")
ab_figobs <- list(
  "treg_ab"=list(
    p_cond_treg, p_ab_fbc_treg, p_ab_test_treg, p_ab_prop_treg
  ),
  "cd4_ab"=list(
    p_cond_cd4, p_ab_fbc_cd4, p_ab_test_cd4, p_ab_prop_cd4
  )
)
# Extremely long running now that we are showing the dotplot.
qs::qsave(ab_figobs, fs::path(outdir, "ab_test_fig.qs"))
