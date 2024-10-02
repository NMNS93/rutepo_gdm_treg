# This script cleans and exports the RUTEPO metadata for publication.
library(readr)
library(dplyr)

## If any fields are updated, additional text files will be required.
meta_txt <- "data/metadata/240813_metadata_cln.csv"
meta <- readr::read_csv(meta_txt)
exp <- meta %>% dplyr::select(
    patient_id, genomic_run, hashtag, weight, height, bmi,
    ethnicity = race, conception,
    has_gdm = gdm_now, gestational_age_weeks = ga_w_out, medication
)
readr::write_csv(exp, "data/export/240813_panicos_to_annotate.csv")
