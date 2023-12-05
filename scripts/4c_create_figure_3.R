source('code/plots.R')

# Setup
library(tidyHeatmap)
library(tibble)
library(dplyr)
library(ggplot2)
figdir = fs::dir_create("gdm/figures_v3/main/")
outdir = 'gdm/results_v3'
# Load DE objects
de_ = qs::qread( fs::path(outdir, "we_treg_cln.qs"))
de_treg = de_$wta; de_cd4 = de_$wca; rm(de_)
de_treg = dplyr::filter(de_treg, !(seurat_clusters %in% c(5,6))) # Remove low cell clusters
# Genes sig DE in any cluster

# Load GSEA objects
tga = qs::qread(fs::path(outdir, "gsea_treg.qs"))
tga = dplyr::filter(tga, !(seurat_clusters %in% c(5,6))) # Remove low cell clusters

# Plot DE dots
de_treg$is_sig = de_treg$is_signif
# Set genes to those significant in eatleast one cluster
tgenes =  de_treg %>% dplyr::filter(is_signif==T) %>% pull(gene) %>% unique()
we_treg = dplyr::filter(de_treg, gene %in% tgenes)
# Set levels of gene so that exclusive genes are ranked
glevels = levels((we_treg %>% group_by(cluster) %>%
  mutate(gene = factor(gene, levels = gene[is_signif])) %>%
  ungroup())$gene)
we_treg$gene = forcats::fct_relevel(we_treg$gene, glevels)
p_de = de_dot_plot(
  we_treg, 
  col.scale=c(gdm_pal[["ogreen"]], "lightgray", gdm_pal[["orange"]])
)

# Plot hallmark GSEA
tgap = tibble(tga) %>%
    dplyr::mutate(
      sig_stars = dplyr::case_when(
          padj <= 0.001 ~ "***",
          padj <= 0.01  ~ "**",
          padj <= 0.05  ~ "*",
          padj <= 0.1 ~ "-",
          TRUE               ~ ""
          ),
      n_edge = sapply(leadingEdge, length),
      pathway_cln = 
        stringr::str_replace_all(
          stringr::str_remove(pathway, "HALLMARK_"),
          "_", " "
      )
)

p_NES = heatmap(
  tibble(tgap),
  cluster, pathway_cln, NES,
  scale = "row",
  palette_value= gdm_pal[c("ogreen", "white", "orange")],
  rect_gp = grid::gpar(col = "#161616", lwd = 0.5)
) %>%
  layer_text(.value=sig_stars, .size=10) %>%
  wrap_heatmap(
    padding = grid::unit(c(0, -80, 0, -100), "points" ), # Aesthetics w/ wrap heatmap b, l, t, r
    clip = FALSE
  )

qs::qsave(list(p_de, p_NES), fs::path(outdir, 'treg_fig3_obs.qs'))

# Create fig 3.
treg_f3 <- patchwork::wrap_plots(
  p_de, p_NES,
  ncol=1,
  heights=c(1,.75)
) + plot_annotation(
  tag_levels=list(c("A", "B"))
)

ggsave(fs::path(figdir, "Figure3.png"), treg_f3, width=8, height=12)



  # LEGACY ----

# # B) Humanin+ genes per-patient Expression and other genes common among clusters
# multi_q_order = dplyr::select(treg@meta.data, multi_q, cond_cln) %>% arrange(cond_cln) %>% pull(multi_q) %>% unique()
# treg$multi_q = factor(treg$multi_q, levels=multi_q_order)
# genes_most_clusters = we_treg_dots %>% dplyr::filter(is_sig==TRUE) %>% group_by(gene) %>% dplyr::count(is_sig) %>% dplyr::filter(n>=5) %>% pull(gene) 
# genes_most_clusters = c('MTRNR2L8', genes_most_clusters)
# genes_most_clusters = genes_most_clusters[order(genes_most_clusters)]
# genes_int = c("MT-CO1", "MT-CO3", "MT-CYB", "MTRNR2L12", "MTRNR2L8")
# DefaultAssay(treg) = "SCT"
# p_common_clust_genes = VlnPlot(treg, genes_int, group.by="multi_q", pt.size=0, stack=T, flip=T, split.by="cond_cln") +
#   scale_fill_manual(values=c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) + theme(legend.position="bottom")
# ggsave(fs::path(figdir, "de_treg_common_viol.png"), p_common_clust_genes, width=5, height=5)
# DefaultAssay(treg) = "SCT"
# # Cache object for use alongside GSEA plot
# de_plots = list(
#   p_treg_dots,
#   p_common_clust_genes
# )
# qs::qsave(de_plots, fs::path(outdir, "de_expression_cache.qs"))

# ### Plot violins for ME of interest
# p_me_vln <- VlnPlot(
#   conn,
#   features = 'Naive-Activation-21',
#   group.by = 'cluster',
#   split.by="cond_cln",
#   pt.size = 0 # don't show actual data points
# ) + geom_boxplot() + scale_fill_manual(values=as.character(gdm_pal[c("orange", "white")]))
# # ============
# p_cd4_dots <- de_dot_plot(
#   we_cd4 %>% get_sig_de_genes, 
#   col.scale=c(gdm_pal[["darkgrey"]], "lightgray", gdm_pal[["pink"]])) +
#   labs(title="CD4")
# ggsave(fs::path(figdir, "de_cd4_dots.png"), p_cd4_dots, width=5, height=12)
# 
# ## C) Stacked violin chart showing significant genes for each cluster
# p_treg_viol_files <- de_viol_plots(
#   we_treg_sig, 
#   order_patients_by_level(treg), "treg",
#   c("lightgray", gdm_pal[["orange"]]), figdir
#   )
# p_cd4_viol_files <- de_viol_plots(
#   we_cd4_sig,
#   order_patients_by_level(cd4), "cd4",
#   c("lightgray", gdm_pal[["pink"]]), figdir)
# 
# # A) Volcanos within clusters
# p_treg_volc <- volcano_within_clusters(we_treg, "treg")
# ggsave(fs::path(figdir, "de_treg_volcano.png"), p_treg_volc, width=15, height=15)
# 
# p_cd4_volc <- volcano_within_clusters(we_cd4, "cd4")
# ggsave(fs::path(figdir, "de_cd4_volcano.png"), p_cd4_volc, width=15, height=15)


# Hallmarks


# HALLMARKS PLOT ----
# Which cluster has the highest enrichment in the most significantly enriched gene set?
# 
# # Create a hallmarks plot for the clusters of interest
# # We focus on Naive-Act-2 
# tg_clusters <- c("Naive-Act-2")
# tgx <- dplyr::filter(tga, cluster %in% tg_clusters & stringr::str_detect(pathway, "HALLMARK")) %>%
#   dplyr::mutate(pathway=stringr::str_remove(pathway, "HALLMARK_")) %>%
#   dplyr::mutate(pathway=stringr::str_replace_all(pathway, "_", " ")) %>% 
#   dplyr::filter(abs(NES)>1)
# p_nes_tgx <- plot_nes_wrap(tgx$cluster, tgx, 
#                            "Hallmark gene set enrichment in\nGDM Naive-Act-2 Tregs",
#                            gdm_pal[["orange"]]) + 
#   theme(legend.position= c(0.85, 0.1), legend.background = element_rect(fill = "white", color = "black", size = 0.5)
#   )
# qs::qsave(p_nes_tgx, fs::path(outdir, "treg_nes.qs"))
# ggsave(fs::path(figdir, "treg_nes.png"), p_nes_tgx, width=7, height=10)
# 
# # GSEA Enrichment plot -----
# 
# # Create an enrichment score plot for the top 3 most significant genes in each cluster
# hpathway = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
# gsea_res <- dplyr::filter(tga, seurat_clusters==6 & pathway == hpathway)
# we_nact <- dplyr::filter(we_treg, seurat_clusters==6 & condition=="CASE")
# p_nfkb_enr = plot_enrich(we_nact, gmts, hpathway, F) + labs(title="TNFA SIGNALING VIA NFKB in\nGDM Naive-Act-2 Tregs")
# qs::qsave(p_nfkb_enr, fs::path(outdir, "treg_nfkb_enr.qs"))
# ggsave(fs::path(figdir, "Treg_Act-2_NFKB_ENRplot.png"), p_nfkb_enr, width=5, height=3)
