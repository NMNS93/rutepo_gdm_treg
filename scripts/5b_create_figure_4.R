# Create Figure 4 plots 

# Setup
source('code/plots.R')
source('code/glm_meta_model.R')
library(data.table)
library(ggplot2)
library(ggupset)
outdir = 'gdm/results_v3/'
figdir = fs::dir_create('gdm/figures_v3/covariates')
maindir = fs::dir_create('gdm/figures_v3/main')

# A: Metadata
meta <- fread(fs::path(outdir, "meta_data.csv"))
p_meta <- plot_gdm_meta(meta)
odds_data = data.table::fread(fs::path(outdir, 'meta_odds.csv'))
p_odds = plot_odds_rat(odds_data %>% dplyr::filter(!stringr::str_detect(variable, "Intercept")))
p_metas <- patchwork::wrap_plots(c(p_meta$indiv, list(p_odds)), nrow=2)
ggsave(fs::path(figdir, "metadata_all.png"), p_metas, width=8, height=8)

# B: Top single markers from all datasets
stirms_smod <- qs::qread(fs::path(outdir, "stirms_smod.qs"))
wang_smod = qs::qread(fs::path(outdir, "wang_smod.qs"))
cd4_smod <- qs::qread(fs::path(outdir, "cd4_smod.qs"))
treg_smod <- qs::qread(fs::path(outdir, "treg_smod.qs"))


# Show markers
smod_all_tops = qs::qread(fs::path(outdir, "smod_all_tops.qs"))

roc_data_af = smod_simple_format(
  smod_all_tops[[2]], smod_all_tops[[1]], smod_all_tops[[3]], smod_all_tops[[4]]
)

roc_data_ss2 <- smod_to_roc_ss(
  smod_all_tops[[2]], smod_all_tops[[1]], smod_all_tops[[3]], smod_all_tops[[4]]
)

prd_f1_plots = plot_roc_f1_single(roc_data_af)

# Show select genes that are in all datasets
roc_data_ss2 <- qs::qread(fs::path(outdir, "roc_data_ss2_cache.qs"))
rdsel = roc_data_ss2 %>% dplyr::select(dataset, model_name) %>% unique() %>% 
  dplyr::mutate(value=1) %>%
  tidyr::pivot_wider(id_cols="model_name", names_from="dataset", values_from='value')

drop_mods = rdsel$model_name[ (rdsel %>% is.na() %>% rowSums() ) > 0]

top_rgenes_bulk = get_top_rgenes_bulk(roc_data_ss2) %>% dplyr::filter(!model_name %in% drop_mods)
prd_roc_plot = plot_roc_data_single(roc_data_ss2 ,
                                    rgenes=top_rgenes_bulk$model_name[1:6], ncol=2)


ggsave(fs::path(figdir, "roc_single.png"), prd_roc_plot,  width=4, height=8)
ggsave(fs::path(figdir, "auc_f1_bulk.png"), prd_f1_plots[[1]], width=8, height=8)
ggsave(fs::path(figdir, "auc_f1_pseu.png"), prd_f1_plots[[2]], width=8, height=8)

# CD4/Treg top 5
## Load Treg top 10
tr5 = qs::qread(fs::path(outdir, "treg_tops5.qs"))$tops
cup = tr5 %>% dplyr::filter(all_y_f1>.5) %>% dplyr::top_n(10, all_y_auc) 
dcup = tops5_to_upset(cup)

## Add CD4 scores at top 10
cd5all = qs::qread( fs::path(outdir, "cd4_mod5.qs")) %>% #Long-running
  dplyr::filter(model_name %in% dcup$model_name)  %>% 
  dplyr::select(model_name, all_cd4_auc=all_y_auc) %>% unique()
dcupa = dcup %>% dplyr::left_join(cd5all, by=c("model_name")) %>%
  tidyr::pivot_longer(c(all_y_auc, all_cd4_auc))

## Gnenerate plot
p_cd4_upset = ggplot(dcupa, aes(x=genes, y=value)) + geom_col(
    aes(fill=name), position=position_dodge()) + scale_x_upset() + theme_bw() + labs(x="", y="CD4+ T pseudobulk\nAUC (LOO-CV)") + scale_y_continuous(n.breaks=10)
p_cd4_upset
ggsave(fs::path(figdir, "cd4_pseu_upset.png"), p_cd4_upset, width=5, height=8)
# 

# RFE boxplots::::
rfe_res = qs::qread(fs::path(outdir, "rfe_singles.qs"))
wang_rfe = rfe_res[[1]]
stirm_rfe = rfe_res[[2]]

wang_rfe_top = wang_rfe[[1]]$optVariables
stirm_rfe_top = stirm_rfe[[1]]$optVariables
wangl = load_wang_bulk("gdm/data/GSE154414_Expression_Gene.csv")
stirm <- load_stirm_bulk("data/mrna/GSE92772_expression_mrna.tsv", "data/mrna/GSE92772_series_matrix.txt")

p_wang_rfe = ggplot(
  dplyr::filter(wangl, gene %in% wang_rfe_top),
  aes(x=gene, y=log10(expression+1), fill=status_c)
) + geom_boxplot(outlier.shape=NA) +
  labs(fill="") + scale_fill_manual(values=c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
  theme_bw() + theme(legend.position="bottom") +
  labs(title="Wang et al. 2021.\nTop CD4+ Markers by RFE") +
  geom_jitter(position=position_jitterdodge(dodge.width=.75, jitter.width=0.1, jitter.height=0.1))

p_stirm_rfe = ggplot(
  dplyr::filter(stirm, gene %in% stirm_rfe_top),
  aes(x=gene, y=expr, fill=status)
) + geom_boxplot(outlier.shape=NA) +
  labs(fill="") + scale_fill_manual(values=c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
  theme_bw() + theme(legend.position="bottom") +
  labs(title="stirm et al. 2018.\nTop CD4+ Markers by RFE") +
  geom_jitter(position=position_jitterdodge(dodge.width=.75, jitter.width=0.1, jitter.height=0.1))

p_rfe = patchwork::wrap_plots(p_wang_rfe, p_stirm_rfe, nrow=1)
ggsave(fs::path(figdir, 'plot_rfe_tops.png'), p_rfe, width=10, height=10)

# ## Figure 8 : DimReg on RFE ----

wdimreg = dplyr::filter(wangl, gene %in% wang_rfe_top) %>% tidyr::pivot_wider(names_from = "gene", values_from="expression")

p_wang_me = ggplot(wdimreg, aes(x=RPL27A, y=TXNIP)) + geom_point(aes(fill=status_c), shape=21, size=7) +
  scale_fill_manual(values=c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
  theme_bw() + theme(legend.position="none") + labs(title="Wang et al., 2021.\nTop Maker Expression")

wang_wide = qs::qread(fs::path(outdir, 'wang_wide.qs'))
l_p_wang = plt_umap_dimreg(wang_wide, wang_rfe_top, 2)
p_wang_umap = l_p_wang[[2]]
p_wang_umap

stw = qs::qread(fs::path(outdir, 'strims_wide.qs'))
l_p_stirm = plt_umap_dimreg(stw, stirm_rfe_top, 5)
p_stirm_umap = l_p_stirm[[2]]
p_stirm_umap
p_stirm_umap = p_stirm_umap + labs(title="Stirm et al., 2018.\nTop CD4+ Marker UMAP") + theme_bw() + theme(legend.position="bottom")

tstt= plt_umap_dimreg(stw,
                      wang_rfe_top %>% purrr::discard(~grepl('MT-ATP8', .x))
                      , 3)

tstt[[1]]
tstt[[2]]
tstt2 = plt_umap_dimreg(
  wang_wide,
  stirm_rfe_top,
  2
)
tstt2[[1]]
tstt2[[2]]

stirm_rfe_top
p_dimreg = patchwork::wrap_plots( p_wang_umap, p_stirm_umap , nrow=1 )

# Seems wang genes define a subset in STirm but Stirm genes do not tranlsate ot wang. 


# Figures ------------
# 
# Cache figure/plot data in case require replotting
figcache = list(
  p_metas,
  p_cd4_upset,
  prd_roc_plot,
  prd_f1_plots,
  p_rfe,
  p_dimreg
)
qs::qsave(figcache, fs::path(outdir, "gdm_covar_figcache.qs"))

# # Create Figure panel ----

# Wrap_elements must be used to be able to label with patchwork
# and to ensure the upset plot x axis does not override
p_panel2 = wrap_elements(full = p_metas[[1]]) +
  wrap_elements(full = p_metas[[4]]) +
  wrap_elements(full = prd_f1_plots[[2]]) +
  wrap_elements(full = prd_f1_plots[[1]]) +
  wrap_elements(full = prd_roc_plot) +
  plot_layout(design = "
              AABBEEE
              CCDDEEE
              CCDDEEE
              ") +
  plot_annotation(tag_levels=c("A"))

  
# TODO: Upset plot as supplementary figure with RFE  

# Note: Cowplot required because upsetplot errors
ggsave(fs::path(maindir, "Figure4.png"), p_panel2, width=14, height=10)

# # Attempt to use plot_grid so that we can align
