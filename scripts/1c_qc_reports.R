# Create preprocessing quality control reports

suppressPackageStartupMessages(source("code/preprocessing_qc.R"))
library(parallel)

outdir = fs::dir_create(here::here("gdm/results_v3"))
figdir = fs::dir_create(here::here("gdm/figures_v3/qc"))
cache_dir = fs::dir_create(fs::path(outdir, "cache"))
sobjs=qs::qread(fs::path(cache_dir, "pre_qc.qs"))
hk_ref = here::here("data/seurat/ref/Housekeeping_GenesHuman.csv")

sobjs = mclapply(sobjs, qc_annotate_tenx, hk_ref=hk_ref, mc.cores=6)

# Plot HK and Mito
hmqc = Reduce(
  rbind,
  lapply(sobjs, function(x) dplyr::select(x@meta.data, percent_mito, percent_hk))
)

p_hmqc = ggplot(hmqc, aes(x=percent_hk, y=percent_mito)) + geom_bin2d() +
  scale_fill_viridis_c() + 
  scale_y_continuous(limits=c(0,25)) +
  scale_x_continuous(limits=c(0,100)) +
  geom_vline(xintercept = 20, linetype="dashed", color="red") +
  geom_hline(yintercept = 5, linetype="dashed", color="red") +
  theme_bw() +
  labs(x="houskeeping genes expressed (%)",
       y="mitochondrial genes expressed (%)")
ggsave(fs::path(figdir, "hmqc.png"), p_hmqc, width=5, height=4)

