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

