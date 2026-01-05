# Export additional data

# Main : Export MTX counts ----

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

# Main : Export Metadata ----

# This script cleans and exports the RUTEPO metadata for publication.
library(readr)
library(dplyr)

## If any fields are updated, additional text files will be required.
meta_txt <- "data/metadata/240813_metadata_cln.csv"
meta <- readr::read_csv(meta_txt)
exp <- meta %>% dplyr::select(
  patient_id, genomic_run, hashtag, weight, height, bmi,
  ethnicity = race, conception,
  gestational_age_weeks = ga_w_out, medication,
)

# Filter to samples in the seurat data object
seu <- qs::qread(fs::path("gdm/results_v3/", "gdm_subtyped.qs"))$treg@meta.data
seu <- dplyr::select(seu, patient_id, has_gdm = cond_cln) %>% unique()
exp2 <- merge(seu, exp)

readr::write_csv(exp, "data/export/251208_panicos_to_annotate.csv")
