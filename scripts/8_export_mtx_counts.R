# Write treg and cd4+ data to format for GEO upload
library(Seurat)
library(fs)
library(glue)
library(readr)
source("code/export_as_tenx.R")

# Load data
gdm <- qs::qread(fs::path("gdm/results_v3/", "gdm_subtyped.qs"))
treg <- gdm$treg
cd4 <- gdm$cd4

# Create mtx file
export_as_tenx(treg, "data/export/treg/")
export_meta(treg, "data/export/treg/")
export_as_tenx(cd4, "data/export/cd4/")
export_meta(cd4, "data/export/cd4")

# Make export directory for GEO
create_geo_directory()
