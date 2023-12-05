# Create figures 1 and 2 from existing data objects. Requires running scripts 1_*->3_*

source("code/plots.R")
source("code/subtyping_de.R")
outdir="gdm/results_v3"
datadir="gdm/data"
figdir =fs::dir_create("gdm/figures_v3/main")
gdm = qs::qread(fs::path(outdir, "gdm_subtyped.qs"))
cd4 = gdm$cd4; treg = gdm$treg; rm(gdm)
saverages = qs::qread(fs::path(outdir, "gdm_average_expr.qs"))
treg_avc = saverages$treg_avc; cd4_avc = saverages$cd4_avc; rm(saverages)
treg_top = readr::read_csv(fs::path(outdir, "top_subtyping_de_markers.csv")) %>%
  dplyr::filter(label=="treg")
cmarkers = readr::read_csv(fs::path(datadir, "paper_markers.csv"))

# Subtype -----

# A
Idents(treg) = "cluster"
treg_tsne = show_tsneplot(treg, pal=treg_dimreg_pal) + labs(title="Treg")

# B
treg_highlight = cmarkers$show[!is.na(cmarkers$show)]
p_trhigh = FeaturePlot(treg, features=treg_highlight, 
                       cols=c(alpha("lightgray",.15), gdm_pal[["orange"]]), reduction="tsne",ncol=3) & 
          theme(legend.position='none', axis.text.x=element_blank(), axis.text.y=element_blank()) 
p_trhigh

# C
top_n = de_top_n_genes(treg_top, 5)
top_papers = cmarkers$genes
genelist = unique(c(top_n, top_papers))

# Reorder gene list by topN in each subtype
dataset = treg_avc@assays$SCT@scale.data 
# Order subtypes for genelist ord by heatmap (based on clust_group)
subtype_order = treg_avc@meta.data %>% select(cluster, clust_group) %>% unique() %>% dplyr::arrange(clust_group) %>% pull(cluster)
genelist_ord <- reorder_by_expr(genelist, treg_avc@assays$SCT@scale.data , (subtype_order), 10)

# Get gene list consisting of top de genes and paper_marker genes.
treg_hm = subtype_heatmap(treg_avc, genelist_ord, is.treg=T)

# treg_f1 <- qs::qread(fs::path(outdir, "pdata_treg_f1.qs"))
# Figure 1
treg_f1 <- patchwork::wrap_plots(
  treg_tsne, p_trhigh, treg_hm + coord_flip(),
  nrow=2,
  heights=c(1,1,0.5),
  design="AAABBB
          AAABBB
          CCCCCC
          CCCCCC"
) + plot_annotation(
  tag_levels=list(c("A", "B", rep("", 8), "C"))
)

ggsave(fs::path(figdir, "Figure1.png"), treg_f1, width=12, height=10)

# Figure 2 -----

# Load abundance figure objects
abfigobs = qs::qread("gdm/results_v3/ab_test_fig.qs")
abtreg = abfigobs$treg_ab

# Set treg pal so it matches figure 1
f1pal_order = c('orange', 'darkred', 'blue', 'ogreen', 'darkgrey', 'purple', 'pink')
f2pal = gdm_dimreg_pal[ match( f1pal_order , names_gdm_dimreg_pal )]

abtreg_2 = abtreg[[2]] + scale_fill_manual(values=f2pal)

treg_f2 <- patchwork::wrap_plots(
  abtreg_2, abtreg[[4]], abtreg[[3]],
  nrow=1,
  widths=c(0.2, 1, 0.4)
) +
  plot_annotation(tag_levels="A") +
  plot_layout(guides='collect') & theme(legend.position='bottom')
ggsave(fs::path(figdir, "Figure2.png"), treg_f2, width=12, height=8)

# 
# # CD4 ----
# # Assignment of T cells to Niave/Effector/Memory subtypes using classical cell surface markers
# 
# # Marker sets, largely derived from Stroukov (Cristiano Paper) and literature
# naivem = c("CCR7", "CD28", "FAS", "SELL") # CD45RA+ (PTPRC alt splicing) --> CCR7+ CD28+ CD95- (FAS)
# tfh = c("CXCR5", "PTPRC")
# thelpers = c("CCR4", "CCR6", "CCR10", "CXCR3") # CXCR3=Th1; CCR6=Th17
# ttreg = c("PTPRC", "FOXP3", "IL2RA", "CTLA4")
# trm2_teff2 = c("CXCR6", "ITGA1", "PRF1", "NKG7") # First 2 are terminal memory, last two are effector memory
# cytotoxic = c("GZMA", "GZMK", "GZMM", "CST7") # GZMK also referred to as effector memory cells
# 
# # Naive subsets: CCR7 is a marker of Naive cells. We observe a set of memory/effector subtypes. FAS negativity further highlights naive cell clusters. CD28 is not very specific. 
# ## Memory or Effector = cluster 2 (c2) and cluster 6 (c6)
# ## Intermediate/Mixed = cluster 4 (c4)
# ## All other clusters Naive
# p_cd4_naive <- subassign_plot(cd4, naivem)
# ggsave(fs::path(figdir, "assign_cd4_naive.png"), p_cd4_naive, width=12, height=8)
# 
# # T follicular Helpers: CXCR5 is a known unique marker, observed in cluster 3. Note that increasing clustering resolution doesn't reveal this cluster as separate.
# ## Mixed = Cluster 3
# p_cd4_tfh <- subassign_plot(cd4, tfh)
# ggsave(fs::path(figdir, "assign_cd4_tfh.png"), p_cd4_tfh, width=12, height=8)
# 
# # T Helpers: Positive for a combination of CCR4 CCR6 CCR10 CXCR3
# # It is clear Cluster 2 carries majority of CCR6 Th cells. The rest are spread among clusters 2 and 3, but not yet clear if these clusters are mix of helpers and effectors
# p_cd4_thelpers <- subassign_plot(cd4, thelpers)
# ggsave(fs::path(figdir, "assign_cd4_thelpers.png"), p_cd4_thelpers, width=12, height=8)
# 
# # Tregs within CD4: We expect a subset of cells to be Tregs, so we search for CD45+ FOXP3+ populations. Although many cells express FOXP3, there is a concentratino in cluster 9.
# ## FOXP3+ = Cluster 9
# p_cd4_tregs = subassign_plot(cd4, ttreg)
# p_cd4_tregs_foxp3 = subassign_thresh(cd4, "FOXP3", .2, ttreg)
# ggsave(fs::path(figdir, "assign_cd4_tregs.png"), p_cd4_tregs, width=12, height=8)
# ggsave(fs::path(figdir, "assign_cd4_tregs_foxp3.png"), p_cd4_tregs_foxp3, width=12, height=8)
# 
# # Memory subtypes: NKG7 shows cluster 6 to be effector memory. This leaves the other CCR7 negative cluster (2) to be memory T cells. 
# ## Effector = Cluster 6
# ## Memory = Cluster 2
# p_cd4_trmems = subassign_plot(cd4, trm2_teff2)
# ggsave(fs::path(figdir, "assign_cd4_trmems.png"), p_cd4_trmems, width=12, height=8)
# p_cd4_trmems_cxcr6 = subassign_thresh(cd4, "CXCR6", .2, trm2_teff2)
# ggsave(fs::path(figdir, "assign_cd4_trmems_cxcr6.png"), p_cd4_trmems_cxcr6, width=12, height=8)
# p_cd4_eff = subassign_plot(cd4, cytotoxic)
# ggsave(fs::path(figdir, "assign_cd4_eff.png"), p_cd4_eff, width=12, height=8)
# 
# ## Information:
# ## Remaining clusters were assigned based on differential expression results.
# ## Cluster 5: C5 is CISH and related to TCR stimulation. This cluster is driven by cells from 1 patient, suggesting immune response at the time of sampling. CISH in T cells: https://assets.researchsquare.com/files/rs-21547/v1/c49a2593-f77d-4582-b11c-cf3ff0cf09b8.pdf?c=1631833841
# ## Cluster 8: MALAT1 promotes terminal effector memory differentiation (https://rupress.org/jem/article/219/6/e20211756/213232/The-long-noncoding-RNA-Malat1-regulates-CD8-T-cell) and EEF1A1 under-expression is associatedw ith cell death (https://www.nature.com/articles/cdd2008136).
# 
# # Treg ----
# 
# # Marker sets (largely from Souztka and literature)
# tr_markers <- c("FOXP3", "IL2RA", "IL7R", "PTPRC") # Stoukas (not FP3) CD25, CD127, CD45RA (splice var).
# tr_thlike <- c("CCR4", "CCR6", "CCR10", "CXCR3") # Same as thelpers
# tr_tfr <- c("CXCR5")
# tr_canonical_sz <- c("FOXP3", "CTLA4", "IRF4", "TNFRSF4") #szabo paper
# tr_prebranch <- c("TCF7", "EEF1B2", "C1orf162", "SNHG7")
# tr_path1e <- c("PTPRC", "DDX17", "MALAT1", "PDE3B")
# tr_path2e <- c("HLA-DRA", "LGALS1", "LGALS3", "CD74")
# 
# # All clusters express FOXP3, except cluster 8
# p_tr_markers = subassign_plot(treg, tr_markers)
# ggsave(fs::path(figdir, "assign_treg_markers.png"), p_tr_markers, width=12, height=8)
# # However, clusters 3, 4 and 9 represent mature (effector) Treg populations with the highest FOXP3, CD25 (prev plot) and CTLA4 expression
# # Interestingly, the other populations do not have low CD127 (IL7R) expression
# p_tr_canonical_sz <- subassign_plot(treg, tr_canonical_sz)
# ggsave(fs::path(figdir, "assign_treg_canonical_sz.png"), p_tr_canonical_sz, width = 12, height = 8)
# 
# # Cluster 7 clear enrichment of CCR6 Th like. Other markers are widespread among the effector path.
# p_tr_thlike <- subassign_plot(treg, tr_thlike)
# ggsave(fs::path(figdir, "assign_treg_thlike.png"), p_tr_thlike, width = 12, height = 8)
# 
# # No cluster with TFR enrichment
# p_tr_tfr <- subassign_plot(treg, tr_tfr)
# ggsave(fs::path(figdir, "assign_treg_tfr.png"), p_tr_tfr, width = 12, height = 8)
# 
# # Pre-branch was identified by TCF7, suggesting 0,1,2,5,6,7 and 10 are precursor states
# p_tr_prebranch <- subassign_plot(treg, tr_prebranch)
# ggsave(fs::path(figdir, "assign_treg_prebranch.png"), p_tr_prebranch, width = 12, height = 8)
# 
# # Cluster 11 showed highest signal for effector Treg path I
# p_tr_path1e <- subassign_plot(treg, tr_path1e)
# ggsave(fs::path(figdir, "assign_treg_path1e.png"), p_tr_path1e, width = 12, height = 8)
# 
# # Cluster 4 showed the highest signal for effector Treg path II
# # Cluster 8 and 4 share LGALS1 expression, but not CD74 or HLA-DRA
# p_tr_path2e <- subassign_plot(treg, tr_path2e)
# ggsave(fs::path(figdir, "assign_treg_path2e.png"), p_tr_path2e, width = 12, height = 8)
# 
# # Assignments summary, based on incorporating DEGs. (Turn this into a wedding notebook)
# # Cluster 0 : Naive-1. DUSP1 under-expression. Low markers of activation like JUN: Naive-Rest1
# # Cluster 1 : Naive-Activation-2. JUN over-expression. Markers of activation. Naive-Act
# # Cluster 2 : Naive-2. Low KLF6 and markers of effector populations such as S100A4
# # Cluster 3 : Effector-Pre. Intermediate state with no canonical marker. Relatively high IL32.
# # Cluster 4 : Effector (Path II). High S100A4 and IL32 expression and CD74 (marker)
# # Cluster 5 : Niave-3. Low expression of activation signal (FOS, etc)
# # Cluster 6 : Naive-4-Activation-1. High expression of activation markers including EGR1, FOS, IER2
# # Cluster 7 : Thelper-Like-1. Expresses CCR6. Also KLRB1. Active genes high with this and low LEF1. Previously described as cytokine producing.
# # Cluster 8 : Effector. ANXA1. Similar cytokine producing population.
# # Cluster 9 : Effector. HLA-A higha nd low EEF1B2 and ARPC1B highest.
# subassign_thresh(treg, "ARPC1B", .2, c("ARPC1B"))
# # Cluster 10 : Naive-5. MTRNR2L12, MTRNR2L8, AL138963.4 expressing. Also carry TCF7 suggesting precursor Treg state. 
# # Cluster 11 : Effector (PathI). High MALAT1, also MTRNR2L12. Low EEF1B2 as with C9 neighbour, could be death of those cells.
# 
# # Load average
# qq = qs::qread(fs::path(outdir, "gdm_average_expr.qs"))
# treg_avc = qq[[1]]
# cd4_avc = qq[[2]]
# rm(qq)
# 
# # Figure 1 ----
# 
# # Heatmap Plots
# Idents(treg) = "cluster"
# treg_tsne = show_tsneplot(treg) + labs(title="Treg")
# treg_hm = subtype_heatmap(treg_avc, de_top_n_genes(dplyr::filter(top_markers, label=="treg"), 5), is.treg=T)
# treg_highlight = c("FOXP3", "CD74", "LGALS1", "TCF7", "S100A4", "JUN")
# p_trhigh = FeaturePlot(treg, features=treg_highlight, cols=c("lightgray", gdm_pal[["orange"]]), reduction="tsne") + NoLegend()
# treg_f1 = patchwork::wrap_plots(list(treg_tsne, p_trhigh, treg_hm), ncols=3) +
#   plot_annotation(tag_levels=list(c(
#     "A", "B", "", "", "", "", "", "C")))
# ggsave(fs::path(figdir, "all_treg_fig1.png"), treg_f1, width=20, height=8)
# # Save plot object for future re-plotting
# qs::qsave(treg_f1, fs::path(outdir, "pdata_treg_f1.qs"))
# 
# Idents(cd4) = "cluster"
# cd4_tsne = show_tsneplot(cd4) + labs(title="CD4")
# cd4_hm = subtype_heatmap(cd4_avc, de_top_n_genes(dplyr::filter(top_markers, label=="cd4"), 5), is.treg=F)
# cd4_highlight = c("CCR7", "FAS", "FOXP3", "GZMA", "S100A4", "JUN")
# p_cd4high = FeaturePlot(cd4, features=cd4_highlight, cols = c("lightgray", gdm_pal[["orange"]]), reduction="tsne") + NoLegend()
# cd4_f1 = patchwork::wrap_plots(list(cd4_tsne, p_cd4high, cd4_hm), ncols=3) 
# ggsave(fs::path(figdir, "all_cd4_fig1.png"), cd4_f1, width=20, height=8)
# # Save plot object for future re-plotting
# qs::qsave(cd4_f1, fs::path(outdir, "pdata_cd4_f1.qs"))
# 
# 
# # Figure_all
# # subtype_fig = patchwork::wrap_plots(cd4_f1, treg_f1, nrow=2) +
# #   plot_annotation(tag_levels=list(c(
# #     "A", "B", "", "", "", "", "", "C",
# #     "D", "E", "", "", "", "", "", "F"))) &
# #     theme(plot.tag=element_text(size=15, face="bold"))
# # ggsave(fs::path(figdir, "all_subtype_fig1.png"), subtype_fig, width=16, height=16)
# 
# 
# # Add case-control contrast
# cd4_cond = show_condition_tsne(cd4) + theme(legend.position="top", text=element_text(size=5)) + labs(title="")
# ggsave(fs::path(figdir, "all_cd4_inset.png"), cd4_cond, width=5, height=5)
# treg_cond = show_condition_tsne(treg) + theme(legend.position="top", text=element_text(size=5)) + labs(title="")
# ggsave(fs::path(figdir, "all_treg_inset.png"), treg_cond, width=5, height=5)
# 
# # # --------------------------
# # # Run scProportionTest
# # # Test for proportions
# # ## N.B. Cluster 8 has IFNG diffexp which is related to diabetes
# # trprop2 <- run_sc_prop_test(treg, "treg")
# # # res = permutation_test(
# # #       sc_utils(treg), cluster_identity = "cluster",
# # #       sample_1 = "CONTROL", sample_2 = "CASE",
# # #       sample_identity = "control"
# # #     )
# # # perm_plt = permutation_plot(res) + labs(title="treg")
# # ggsave(fs::path(figdir, "permutation_treg1.png"),trprop@plot, width=10, height=5)
# # ggsave(fs::path(figdir, "permutation_treg2.png"),trprop2@plot, width=10, height=5)
# # # qs::qsave(res, fs::path(outdir, "treg_permtest.qs"))
# # 
# # cd4prop <- run_sc_prop_test(cd4, "cd4")
# # ggsave(fs::path(figdir, "permutation_cd4.png"),cd4prop@plot, width=10, height=5)
# # 


# Abundance -----

# Plot figures together
# Create Figure 2 for GDM abundance ----
# 
# panel_a <- patchwork::wrap_plots(
#   ab_figobs$cd4_ab,
#   nrow=1,
#   widths=c(0.7,0.5,1.3)
# ) +
#   plot_annotation(tag_levels=list(c(
#     "A", "B", "C", "D")))
# ggsave(fs::path(figdir, "abundance_f1a.png"), panel_a, width=15, height=8) 
# 
# panel_b <- patchwork::wrap_plots(
#   ab_figobs$treg_ab,
#   nrow=1,
#   widths=c(0.7,0.5,1.3)
# ) +
#   plot_annotation(tag_levels=list(c(
#     "A", "B", "C", "D"))) 
# ggsave(fs::path(figdir, "abundance_f2a.png"), panel_b, width=15, height=8)  
# 
# #panels_fig2 <- patchwork::wrap_plots(panel_a, panel_b, nrow=2, ncol=1)
# #ggsave(fs::path(figdir, "abundance_fig2.png"), panels_fig2, width=15, height=18)
# 
# # Create composite plot for Figure 1
# treg_f1 <- qs::qread(fs::path(outdir, "pdata_treg_f1.qs"))
# panel_treg <- patchwork::wrap_plots(
#   list(treg_f1, panel_b),
#   nrow=2
#   #widths=c(0.7,0.5,1.3)
# ) +
#   plot_annotation(tag_levels=list(c(
#     "A", "B", "", "", "", "", "", "C",
#     "D", "E", "F", "G"))) &
#   theme(plot.tag=element_text(size=25, face="bold"))
# ggsave(fs::path(figdir, "PanelTregF1.png"), panel_treg, width=18, height=15, dpi=500)
# 
# ## TODO: Repeat for CD4
# ## 