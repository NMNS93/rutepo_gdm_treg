# Assign subtypes by updating the labels_v*.csvs prior to running this script. 
# Subtypes informed by figures created in 2_subtyping.R

# Setup
source("code/plots.R")
source("code/subtyping_de.R")
outdir=fs::dir_create("gdm/results_v3/")
figdir=fs::dir_create("gdm/figures_v3/subtyping")
top_markers = readr::read_csv(fs::path(outdir, "top_subtyping_de_markers.csv"))
gdm = qs::qread(fs::path(outdir, "gdm_seurat_sct_markers.qs"))
treg = gdm$treg; cd4=gdm$cd4; rm(gdm)
DefaultAssay(treg) <- "SCT"
DefaultAssay(cd4) <- "SCT"

# Manually updated data files used to assign labels:
cd4_labels <- "gdm/data/cd4_labels_v3.csv"; file.edit(cd4_labels)
treg_labels <- "gdm/data/treg_labels_v3.csv"; file.edit(treg_labels)

# Assign subtypes to objects and save
cd4 <- add_cluster_labels(cd4, cd4_labels)
treg <- add_cluster_labels(treg, treg_labels)
qs::qsave(list("cd4"=cd4, "treg"=treg), fs::path(outdir, "gdm_subtyped.qs"))
gdm = qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 = gdm$cd4
treg = gdm$treg
rm(gdm)

# Save average expression per cluster for each gene (for heatmap)
cd4_avc = get_cluster_average(cd4)
cd4_avc = add_label_data(cd4_avc, cd4_labels)
treg_avc = get_cluster_average(treg)
treg_avc = add_label_data(treg_avc, treg_labels)
qs::qsave(
  list("treg_avc"=treg_avc, "cd4_avc"=cd4_avc),
  fs::path(outdir, "gdm_average_expr.qs")
)
