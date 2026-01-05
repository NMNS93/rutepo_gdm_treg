# Generate visualisations

# Library ----
library(ggplot2)
library(Seurat)
library(glue)
library(patchwork)
library(ggsignif)
library(umap)
library(forcats)
library(ggrepel)
library(dplyr)

# Palettes ----

# Project colour palette
gdm_pal <- c(
  "black" = "#000000",
  "darkgrey" = "#252525", "darkpurp" = "#171723",
  "ogreen" = "#004949", "purple" = "#490092",
  "darkred" = "#920000", "brown" = "#8f4e00",
  "green" = "#22cf22", "white" = "#ffffff",
  "lightgrey" = "#676767", "blue" = "#006ddb", "teal" = "#009999",
  "lightpurple" = "#b66dff",
  "pink" = "#ff6db6", "orange" = "#db6d00", "yellow" = "#ffdf4d"
)

# Project colour palette with alpha
gdm_dimreg_pal <- as.character(alpha(
  gdm_pal[c(
    "purple", "lightpurple", "blue", "darkred",
    "pink", "yellow", "teal", "orange", "ogreen",
    "green", "brown", "darkgrey", "darkpurp"
  )],
  alpha = 0.8
))
names_gdm_dimreg_pal <- c(
  "purple", "lightpurple", "blue", "darkred",
  "pink", "yellow", "teal", "orange", "ogreen",
  "green", "brown", "darkgrey", "darkpurp"
)

# Condition palettes
cond_cln_pal <- as.character(alpha(c("red", "blue"), alpha = 0.3))
treg_pal <- c("purple", "pink", "blue", "darkred", "orange", "ogreen", "darkgrey")
treg_dimreg_pal <- gdm_dimreg_pal[match(treg_pal, names_gdm_dimreg_pal)]

# Dataset palettes
cd4_pal <- c(
  "purple", "lightpurple", "blue", "darkred",
  "pink", "darkgrey", "teal", "orange", "ogreen",
  "green", "brown", "yellow", "darkpurp"
)
cd4_dimreg_pal <- gdm_dimreg_pal[match(cd4_pal, names_gdm_dimreg_pal)]


# Function ----

#' Visualisate clusters
show_cluster <- function(x, cv, tt = "", reduction = "tsne", col = "red") {
  ccells <- WhichCells(x, expression = seurat_clusters == cv)
  DimPlot(x, cols.highlight = col, cells.highlight = ccells, reduction = reduction, sizes.highlight = 0.5) + labs(title = tt) +
    NoLegend()
}

#' Visualise features
show_feature <- function(x, feature, reduction = "tsne") {
  FeaturePlot(x, feature, reduction = reduction)
}

#' Iterate over features
iter_show_feature <- function(x, features, reduction = "tsne") {
  for (i in features) {
    print(show_feature(x, i, reduction))
    invisible(readline(prompt = ""))
  }
}

#' Visualise counts per cluster per patient
show_count_pat_clust <- function(x) {
  ggplot(x@meta.data, aes(x = forcats::fct_reorder(multi_q, control == "CASE"), fill = control)) +
    geom_bar(color = "black", size = 0.5) +
    scale_fill_viridis_d() +
    facet_wrap(~seurat_clusters) +
    theme_classic() +
    theme(axis.text.x = element_blank()) +
    labs(x = "", y = "N cells")
}

#' Visualise features
show_feat <- function(x, feat) {
  FeaturePlot(x, feat, reduction = "tsne") / VlnPlot(x, feat)
}

#' Plot median cells per patient per cluster
plot_med_per_pat_clust <- function(x) {
  ggplot(x, aes(x = seurat_clusters, y = med)) +
    geom_point() +
    geom_errorbar(aes(ymin = med - (iqr / 2), ymax = med + (iqr / 2), width = 0.15)) +
    theme_bw() +
    labs(title = x$label[[1]], y = "Medain number of cells per patient") +
    scale_y_continuous(n.breaks = 20)
}

#' Composite dimensionality reduction plots
show_dimplots <- function(x, y) {
  px <- (UMAPPlot(x, cols = gdm_dimreg_pal) + NoLegend()) | (TSNEPlot(x, cols = gdm_dimreg_pal))
  py <- (UMAPPlot(y, cols = gdm_dimreg_pal) + NoLegend()) | (TSNEPlot(y, cols = gdm_dimreg_pal))
  px / py
}

#' Visualise TSNE plot
show_tsneplot <- function(x, pal = gdm_dimreg_pal) {
  plt <- (TSNEPlot(x, cols = pal)) + NoLegend() # + theme(legend.position="bottom")
  LabelClusters(plt, id = "ident", box = T, color = "white")
}

#' Visualise TSNE per condition
plot_cond_tsne <- function(x, col_gdm) {
  Idents(x) <- "cond_cln"
  TSNEPlot(x, cols = c(col_gdm, gdm_dimreg_pal[12])) + theme(legend.position = "bottom") + labs(title = "")
}

#' Add custom border theme
theme_cust_border <- function() {
  theme(
    panel.border = element_rect(fill = NULL, color = "black"),
    axis.line = element_blank()
  )
}

#' Visualise SCT tsne clusters
ShowSCTClusterDeTop <- function(x, de, reduction = "tsne") {
  DefaultAssay(x) <- "SCT"
  cc <- sort(unique(x$seurat_clusters))

  # Add color based on annotation
  de$color <- ifelse(de$CONTROL_avg_log2FC >= 0, "red", "blue")

  # Rank de top 3 genes based on fold-change. All significant at this point.
  de3 <- dplyr::group_by(de, ident) %>% dplyr::top_n(5, abs(CONTROL_avg_log2FC))

  # Split de by cluster
  desplit <- split(de3, de3$ident)

  # Loop over cluster and generate plot for all three top genes
  deplots <- lapply(
    names(desplit),
    function(cluster) {
      data <- desplit[[cluster]]
      p_highlight <- show_cluster(x, cluster, cluster, reduction = reduction)
      p_tops <- lapply(
        seq(nrow(data)),
        function(i) {
          gene <- data[i, ]$gene
          col <- data[i, ]$color
          FeaturePlot(x, gene, reduction = reduction, cols = c("lightgrey", col))
        }
      )
      patchwork::wrap_plots(c(list(p_highlight), p_tops), ncol = 4)
    }
  )

  # Use patchwork to combine plots into single column plot
  p_topw <- patchwork::wrap_plots(deplots, ncol = 1)
  return(p_topw)
}

#' Show differential expression results as dot plot
de_dot_plot <- function(x, col.scale = c("blue", "grey", "red")) {
  # x = DE results table
  # Filter for significant genes
  x$is_sig <- ifelse(x$is_sig == T & !is.na(x$is_sig), T, F)
  ggplot(x, aes(y = gene, x = cluster, fill = avg_log2FC, size = -log10(p_val_adj))) +
    geom_point(aes(shape = is_sig, alpha = is_sig)) +
    # geom_point(aes(fill=avg_log2FC, size=-log10(p_val_adj)*0.75), shape=21, alpha=0.5) +
    # geom_point(aes(fill=avg_log2FC, size=-log10(p_val_adj)), shape=22, data=dplyr::filter(x, is_sig==T)) +
    theme_bw() +
    labs(x = "", y = "") +
    scale_x_discrete(position = "bottom") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1.05, vjust = 1.05)) +
    scale_fill_gradient2(
      low = col.scale[[1]], mid = col.scale[[2]], high = col.scale[[3]],
      n.breaks = 8
    ) +
    scale_color_manual(
      values = c("FALSE" = "white", "TRUE" = "black"),
      breaks = c(TRUE, FALSE), # Values to include in the legend
      labels = c("q-value ≤ 0.05", "non-significant") # Labels to replace TRUE and FALSE
    ) +
    scale_alpha_discrete(range = c(0.55, 1), labels = c("N.S.", "adjusted p ≤ 0.05")) +
    scale_shape_manual(values = c(21, 22), labels = c("N.S.", "adjusted p ≤ 0.05")) +
    labs(shape = "Significance", size = "-log10\n(adjusted p)", fill = "Average Log2\nFold-change", alpha = "Significance")
}

#' Helper.
subassign_plot <- function(x, features) {
  (FeaturePlot(x, features, col = c("white", "blue"), reduction = "tsne") | TSNEPlot(x, cols = gdm_dimreg_pal)) /
    VlnPlot(x, features, pt.size = 0, cols = gdm_dimreg_pal)
}

#' Helper.
subassign_thresh <- function(x, feature, thresh, features) {
  x <- subset(x, !!sym(feature) >= thresh)
  subassign_plot(x, features) + plot_annotation(title = paste0(feature, ">=", thresh))
}

#' Helper.
subtype_heatmap <- function(avg_per_clust, genelist, is.treg = T) {
  hm <- DoHeatmap(avg_per_clust,
    size = 5, group.bar = F, # group.by="clust_group",
    features = genelist, angle = 0, hjust = 0.3,
    assay = "SCT", draw.lines = F,
    group.colors = gdm_dimreg_pal[c(1, 6, 10)]
  ) + guides(color = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  if (is.treg) {
    hm <- hm +
      scale_fill_gradientn(colors = c(gdm_pal[["teal"]], gdm_pal[["white"]], gdm_pal[["orange"]]))
  } else {
    hm <- hm +
      # Make consistent as colours for CD4 and Treg are not informative
      # scale_fill_gradientn(colors=c(gdm_pal[["darkgrey"]], gdm_pal[["white"]], gdm_pal[["pink"]]))  +
      scale_fill_gradientn(colors = c(gdm_pal[["teal"]], gdm_pal[["white"]], gdm_pal[["orange"]]))
  }
  return(hm)
}

#' Helper.
show_condition_tsne <- function(x) {
  TSNEPlot(x, group.by = "cond_cln", cols = cond_cln_pal)
}

#' Visualise abundance
plot_abundance_fbc <- function(x, tt, pal = gdm_dimreg_pal) {
  ggplot(x@meta.data, aes(x = cond_cln)) +
    geom_bar(aes(fill = cluster), position = "fill", color = "black") +
    scale_fill_manual(values = pal) +
    theme_classic() +
    labs(x = "", y = "Proportion", fill = "", title = tt) +
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
}

#' Visualise propeller results
plot_propeller_res <- function(x, tt, sig_col) {
  ggplot(x, aes(y = cluster, x = PropRatio)) +
    geom_vline(xintercept = c(.75, 1.25), linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 1, color = "black") +
    geom_point(size = 3.5) +
    geom_point(color = sig_col, size = 2, data = x %>% dplyr::filter(abs(FDR) <= 0.1)) +
    theme_bw() +
    guides(color = "none") +
    labs(x = "Proportion Ratio (GDM:Control)", y = "", title = "") +
    scale_x_continuous(n.breaks = 6)
}

#' Visualise patient proportions
plot_pat_proportion <- function(x, col) {
  case_levels <- x@meta.data %>%
    dplyr::select(multi_q, cond_cln) %>%
    dplyr::arrange(cond_cln, multi_q) %>%
    unique() %>%
    dplyr::pull(multi_q)

  pdata <- x@meta.data %>%
    dplyr::group_by(multi_q) %>%
    dplyr::mutate(pat_cells = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(multi_q, cond_cln, cluster, pat_cells) %>%
    dplyr::summarise(clust_cells = n()) %>%
    dplyr::mutate(clust_prop = clust_cells / pat_cells) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(multi_q = fct_relevel(multi_q, case_levels))

  ggplot(pdata, aes(x = multi_q, y = clust_prop, fill = cond_cln)) +
    geom_col() +
    theme_bw() +
    scale_y_continuous(n.breaks = 3) +
    facet_wrap(~ fct_rev(cluster), ncol = 1, scales = "free_y", strip.position = "right") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = c("black", col)) +
    theme(
      strip.background = element_rect(fill = "white", color = NULL),
      legend.position = "top", # strip.text=element_text(size=4),
      axis.text.x = element_blank()
    ) +
    labs(x = "", fill = "", y = "Proportion")
}

#' Plot GDM metadata
plot_gdm_meta <- function(data) {
  p_ethnicity <- ggplot(data, aes(x = ethnicity, fill = status_c)) +
    geom_bar(position = "dodge", color = "black", alpha = 0.5) +
    labs(title = "", x = "Ethnicity", y = "Count") +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    theme_bw() +
    scale_y_continuous(n.breaks = 10) +
    theme(legend.position = "none") +
    labs(fill = "", y = "")

  p_bmi <- ggplot(data, aes(x = status_c, y = bmi, fill = status_c)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(alpha = 0.5, width = 0.13, shape = 21) +
    labs(title = "", x = "", y = "Body mass index (BMI)") +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    scale_y_continuous(limits = c(min(data$bmi) - 5, NA), n.breaks = 8) +
    theme_bw() +
    theme(legend.position = "none")

  p_map <- ggplot(data, aes(x = status_c, y = map, fill = status_c)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(alpha = 0.5, width = 0.13, shape = 21) +
    labs(title = "", x = "", y = "Mean arterial pressure (MAP)") +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    scale_y_continuous(limits = c(min(data$map) - 5, NA), n.breaks = 8) +
    theme_bw() +
    theme(legend.position = "none")

  p_meta <- patchwork::wrap_plots(p_ethnicity, p_bmi, p_map, ncol = 1)
  return(list(plot = p_meta, indiv = list(p_bmi, p_ethnicity, p_map)))
}

#' Plot bulk mRNA boxplot
plot_stirms_box <- function(strims, figdir) {
  bulk_plot <- ggplot(strims, aes(x = status, y = expr)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~gene)
  ggsave(fs::path(figdir, "strims_boxplots.png"), bulk_plot, width = 10, height = 10)
}

#' Bulk mRNA T-test
plot_stirms_gene <- function(strims) {
  ggplot(strims, aes(x = status, y = expr, fill = status)) +
    geom_violin(alpha = 0.8) +
    stat_summary(fun = "median", geom = "point", shape = "-", size = 40, alpha = 0.5) +
    geom_signif(comparisons = list(c("GDM", "NGT"))) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = strims$gene[1]) +
    scale_fill_manual(values = c(gdm_pal[["orange"]], gdm_pal[["lightgrey"]]))
}

#' Helper.
p_wrap_stirms_gene <- function(strims, figdir) {
  strims_genes <- unique(strims$gene)
  strims_gene_plots <- lapply(
    unique(strims$gene),
    function(x) {
      plot_stirms_gene(dplyr::filter(strims, gene == x))
    }
  )
  p_sgp <- patchwork::wrap_plots(strims_gene_plots, ncol = 6)
  ggsave(fs::path(figdir, "strims_violins_pval.png"), p_sgp, width = 20, height = 20)
  return(p_sgp)
}

#' Plot CD4 top genes
plot_cd4_top_genes_box <- function(dbox) {
  dbox %>%
    ggplot(aes(x = gene, y = expr, fill = status_c)) +
    geom_violin() +
    stat_summary(
      fun = median, fun.min = median, fun.max = median,
      color = gdm_pal[["yellow"]], geom = "crossbar", width = 0.3, position = position_dodge(width = 0.9)
    ) +
    theme_bw() +
    labs(x = "", y = "Log10(CPM)") +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    geom_jitter(
      shape = 21,
      position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.05),
      alpha = 0.5
    ) +
    facet_wrap(~dataset, ncol = 1) +
    theme_light()
}

#' Plot CD4 UMAP
plt_umap_dimreg <- function(cd4_dataset, gtop, neighbors = 5, ptitle = "") {
  pca_result <- prcomp(cd4_dataset[, gtop], scale = F)

  # Extract PCA scores
  pca_scores <- as.data.frame(pca_result$x[, 1:2])
  pca_scores$status <- ifelse(cd4_dataset$status == 1, "GDM", "CONTROL")

  # Plot PC1 and PC2 using ggplot
  plt_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = status)) +
    geom_point() +
    labs(title = ptitle, x = "Principal Component 1", y = "Principal Component 2") +
    theme_minimal() +
    scale_color_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]]))

  umap_result <- umap(cd4_dataset[, gtop], n_neighbors = neighbors)

  # Combine UMAP results with status column
  umap_data <- data.frame(cbind(umap_result$layout, data.frame(status = ifelse(cd4_dataset$status == 1, "GDM", "CONTROL"))))

  # Create a scatter plot using ggplot
  plt_umap <- ggplot(umap_data, aes(x = X1, y = X2, fill = status)) +
    geom_point(size = 7, shape = 21, color = "black") +
    labs(title = ptitle, x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
    theme_minimal() +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]]))

  plt_dimreg <- patchwork::wrap_plots(plt_umap, plt_pca, nrow = 1) + plot_layout(guides = "collect")
  plt_dimreg
  return(list(plt_pca, plt_umap))
}

#' Plot ROC data
plot_roc_data_single <- function(roc_data_ss, rgenes, ncol = 2) {
  rss <- dplyr::filter(roc_data_ss, model_name %in% rgenes)
  udatasets <- unique(rss$dataset)
  # Plot ROC curves
  roc_text <- dplyr::select(rss, dataset, model_name, all_y_auc) %>%
    unique() %>%
    dplyr::mutate(label = paste0(dataset, "-AUC:", round(all_y_auc, 2))) %>%
    dplyr::mutate(y = dplyr::case_when(
      dataset == udatasets[4] ~ 0.4,
      dataset == udatasets[1] ~ 0.3,
      dataset == udatasets[2] ~ 0.2,
      dataset == udatasets[3] ~ 0.1
    ), x = 0.7)
  p_strim_roc <- ggplot(rss, aes(m = y_prob, d = y_truth, color = dataset)) +
    geom_roc(n.cuts = 0, alpha = 0.5) +
    geom_text(data = roc_text, mapping = (aes(x = x, y = y, label = label, color = dataset)), inherit.aes = F) +
    facet_wrap(~model_name, ncol = ncol) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_colour_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["teal"]], gdm_pal[["pink"]], gdm_pal[["orange"]])) +
    labs(x = "False Positive Fraction", y = "True Positive Fraction") +
    theme(strip.background = element_rect(fill = "white"))
  return(p_strim_roc)
}

#' Plot ROC data with F1 score
plot_roc_f1_single <- function(x) {
  df <- x %>%
    dplyr::select(-fold_no, -fold_id, -y_truth, -y_prob) %>%
    unique() # %>% dplyr::filter(all_y_auc >= 0.5 & all_y_f1 >=0.5)
  df <- df %>%
    mutate(rank_metric = all_y_f1 + all_y_auc) %>%
    arrange(desc(rank_metric))
  df$is_pseu <- ifelse(
    stringr::str_detect(df$dataname, "scRNA"),
    "Single-cell Pseudobulk",
    "Bulk RNA"
  )
  datasets <- split(df, df$is_pseu)

  plt <- ggplot(datasets[[1]], aes(x = all_y_f1, y = all_y_auc)) +
    geom_line(aes(group = model_name), linetype = "dashed", alpha = 0.6) +
    geom_jitter(aes(color = dataname), size = 3, alpha = .5, shape = 17) +
    geom_label_repel(aes(label = model_name),
      max.overlaps = 50,
      alpha = 0.8, box.padding = 0.25, size = 3.5, force = 10, max.iter = 200
    ) +
    theme_bw() +
    scale_color_manual(values = c(as.character(gdm_pal[c("green", "orange")]))) +
    labs(x = "F1 Score", y = "AUC Value", color = "Bulk RNA") +
    theme(legend.position = c(0.75, 0.15), legend.box.background = element_rect(color = "black", size = 0.5)) +
    scale_x_continuous(limits = c(0.25, 1)) +
    scale_y_continuous(limits = c(0.25, 1))


  plt2 <- ggplot(datasets[[2]], aes(x = all_y_f1, y = all_y_auc)) +
    geom_line(aes(group = model_name), linetype = "dashed", alpha = 0.6) +
    geom_jitter(aes(color = dataname), size = 3, alpha = .5) +
    geom_label_repel(aes(label = model_name),
      max.overlaps = 50,
      alpha = 0.8, box.padding = 0.25, size = 3.5
    ) +
    theme_bw() +
    scale_color_manual(values = c(as.character(gdm_pal[c("darkgrey", "pink")]))) +
    labs(x = "F1 Score", y = "AUC Value", color = "Pseudobulk") +
    theme(legend.position = c(0.75, 0.15), legend.box.background = element_rect(color = "black", size = 0.5)) +
    scale_x_continuous(limits = c(0.25, 1)) +
    scale_y_continuous(limits = c(0.25, 1))

  return(list(plt, plt2))
}

#' Plot rat odds ratio
plot_odds_rat <- function(x) {
  # Create a ggplot for odds ratios and CIs
  ggplot(x, aes(x = odds_ratio, y = clean, xmin = lower_ci, xmax = upper_ci)) +
    geom_point(stat = "identity", position = position_dodge(width = 0.75), size = 3) +
    geom_errorbar(position = position_dodge(width = 0.75), width = 0.2) +
    labs(
      x = "Odds Ratio",
      y = ""
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))
}

#' Top 5 upset plots
tops5_to_upset <- function(x) {
  cup <- x %>% dplyr::select(all_y_auc, model_name)
  ugenes <- unique(unlist((stringr::str_split(cup$model_name, ";"))))
  ugenes <- ugenes[ugenes != "status"]
  ugenes_bool <- sapply(ugenes, function(x) stringr::str_detect(cup$model_name, x))
  cup <- cbind(cup, ugenes_bool)
  cup$status <- NULL
  cupl <- tidyr::pivot_longer(cup, cols = all_of(ugenes), names_to = "gene", values_to = "present")
  dcup <- cupl %>%
    filter(present == TRUE) %>%
    group_by(model_name, all_y_auc) %>%
    summarise(genes = list(gene)) %>%
    dplyr::arrange(-all_y_auc)
  return(dcup)
}
