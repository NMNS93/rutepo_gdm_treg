# Run GSEA on hallmarks for all scRNA clusters with DE
#   Use both p_values and log2FC as ranking

source("code/gsea.R")
source("code/plots.R")
library(qs)
library(fs)
library(dplyr)
library(readr)
library(logger)
library(parallel)
library(gridExtra)

# Load data
outdir = here::here("gdm/results_v3")
figdir = fs::dir_create("gdm/figures_v3/gsea")
gdm = qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 = gdm$cd4; treg = gdm$treg; rm(gdm)

# Load differential expression
we <- qs::qread(fs::path(outdir,"we_treg_cln.qs"))
we_treg = we$wta; we_cd4 = we$wca

# Load hallmark GMT files for initial GSEA
gmts <- fs::dir_ls("data/gsea_gene_sets/", glob="*h.all.v7.5.1.symbols.gmt")

# Run GSEA for each cluster
log_info("Treg gsea")
we_treg_pass <- dplyr::filter(we_treg, pct.1>=0.05 & pct.2 >=0.05)
tga <- gsea_runner_v3(we_treg_pass, gmts)
qs::qsave(tga, fs::path(outdir, "gsea_treg.qs"))

log_info("CD4 gsea")
we_cd4_pass <- dplyr::filter(we_cd4, pct.1>=0.05 & pct.2 >=0.05)
cga <- gsea_runner_v3(we_cd4, gmts)
qs::qsave(cga, fs::path(outdir, "gsea_cd4.qs"))



# Founder gene sets ----

# For the most significant hallmark pathways in clusters of interest, perform GSEA on the founder gene sets
# # TOP Hallmarks: TNFA signalling, UV response, IFN gamma, Hypoxia. GMTs downloaded manually.
# hgmts <- fs::dir_ls("data/gsea_gene_sets/", glob="*HALLMARK*.gmt")
# log_info("Treg HGMT")
# htga <- gsea_runner(we_treg, hgmts, field="avg_log2FC")
# qs::qsave(htga, "temp.qs") #del
# log_info("CD4_HGMT")
# hcga <- gsea_runner(we_cd4, hgmts, field="avg_log2FC")
# qs::qsave(hcga, "temp2.qs") #del
# log_info("Saving HGMT")
# qs::qsave(
#   list("htga"=htga, "hcga"=hcga), # T and C for Treg and CD4 respectively
#   fs::path(outdir, "results_fgsea_founders.qs")
# )
# htgdata = qs::qread(fs::path(outdir, "results_fgsea_founders.qs"))
# 
# 
# # Identify the gene sets contributing to the Act-2 TNFA hallmark signal
# tnf_founder = htga$gmt %>% purrr::keep(~stringr::str_detect(., "TNF")) %>% unique()
# top_foundersx = dplyr::filter(htga, gmt==tnf_founder & seurat_clusters==6 & padj <= 0.01 & size>10 & abs(NES) > 2) %>% arrange(abs(padj))
# top_founders = top_foundersx$pathway[1:5]

## SS-GSEA ----
# ## Run ssgsea using ESCAPE to visualise the NFKB gene set across clusters 
# ### First, select gene sets for each dataset
# gml = gmtPathways(paste0("data/gsea_gene_sets/", tnf_founder))
# gene.sets = gml[top_founders]
# 
# ### Next, run ESCAPE on each
# library(escape) # Long-loading library
# tES <- enrichIt(obj = treg, 
#                 gene.sets = gene.sets, 
#                 groups = 1000, cores = 8, 
#                 min.size = 10)
# qs::qsave(list("tES"=tES), fs::path(outdir, "escape_scGSEA_v2.qs"))
# escp = qs::qread(fs::path(outdir, "escape_scGSEA_v2.qs"))
# 
# ### Add metadata
# tss = AddMetaData(treg, tES)
# 
# ### Run ESCAPE significance testing
# ### These gene sets were already selected using fGSEA significance testing on Naive-Act-2.
# ### ESCAPE T-test should show that these, like the hallmarks, remain most-significant on Naive-Act-2.
# tssl <- tss@meta.data %>% dplyr::select( c( all_of(top_founders), cluster, cond_cln ) )
# tssl = split(tssl, tssl$cluster)
# tssig <- lapply(
#   names(tssl),
#   function(x){
#     cur = tssl[[x]]
#     res = escape::getSignificance(cur, group="cond_cln", fit="T.test")
#     res$cluster = x
#     return(res)
#   }
# )
# tssig = Reduce(rbind, tssig)
# tssig$padj = p.adjust(tssig$p.value, method = "BH")
# tssig$gene_set = rownames(tssig) 
# tssig$gene_set = gsub("\\d+$", "", tssig$gene_set)
# rownames(tssig) = NULL
# readr::write_csv(tssig, fs::path(outdir, "ESCAPE_TopNaive-Act-2_Ttest.csv"))
# 
# # Plot hallmark founders
# thfplot <- p_hallmark_founder_sig(tssig, "Top Hallmark TNFA signalling via NFKB Founder Gene Sets (scGSEA)")
# ggsave(fs::path(figdir, "TNFA_Nact-2_hallmarks.png"), thfplot, width=10, height=5)
# 
# # Next, build a heatmap using the Naive-Act2 cells only
# library(dittoSeq)
# tss_nact2 = subset(tss, seurat_clusters == 6)
# p_treg_esc_hm = build_escape_heatmap(
#   tss_nact2,
#   #fs::path(figdir, "treg_nact2_escape_heatmap.png"),
#   genesets=top_founders,
#   c(gdm_pal["green"], "white", gdm_pal["orange"])
# )
# qs::qsave(p_treg_esc_hm, fs::path(outdir, "p_treg_escape_hm.qs"))
# ggsave(fs::path(figdir, "treg_nact2_escape_heatmap.png"), p_treg_esc_hm, width=10, height=3)

# Future option: p_treg_esc_hm may be recreated with heatmaply directly as ggplot2 object if all else fails

# ### RidgeEnrichment
# for(i in names(htsig[[2]])){
#   print(paste0("Making escape:ridge+split:", i))
#   make_escape_ridge_split(tss, i, c(gdm_pal["orange"], "white", gdm_pal["green"]), "treg", figdir)
# }

## Note: SIRT1 and NFKB are related. Piccaluga AIL is known to be linked with NFKB downreg.