library(ggplot2)
library(Seurat)
library(glue)
library(patchwork)
library(ggsignif)
library(umap)
library(forcats)
library(ggrepel)
library(dplyr)

# GDM Project color palette: https://github.com/filipworksdev/colorblind-palette-16
gdm_pal = c("black"="#000000",
        "darkgrey"="#252525","darkpurp"="#171723",
        "ogreen"="#004949","purple"="#490092",
        "darkred"="#920000","brown"="#8f4e00",
        "green"="#22cf22", "white"="#ffffff", 
        "lightgrey"="#676767", "blue"="#006ddb", "teal"="#009999",
        "lightpurple"="#b66dff", 
        "pink"="#ff6db6", "orange"="#db6d00", "yellow"="#ffdf4d"
)
gdm_dimreg_pal <- as.character(alpha(
  gdm_pal[c("purple", "lightpurple", "blue", "darkred",
            "pink", "yellow", "teal", "orange", "ogreen",
            "green", "brown", "darkgrey", "darkpurp")],
  alpha=0.8
))
names_gdm_dimreg_pal <- c("purple", "lightpurple", "blue", "darkred",
                           "pink", "yellow", "teal", "orange", "ogreen",
                           "green", "brown", "darkgrey", "darkpurp")
cond_cln_pal <- as.character(alpha(c("red","blue"), alpha=0.3))
treg_pal <- c("purple", "pink", "blue", "darkred", "orange", "ogreen", "darkgrey")
treg_dimreg_pal <- gdm_dimreg_pal[ match(treg_pal, names_gdm_dimreg_pal )]

cd4_pal <- c("purple", "lightpurple", "blue", "darkred",
             "pink", "darkgrey", "teal", "orange", "ogreen",
             "green", "brown", "yellow", "darkpurp")
cd4_dimreg_pal <- gdm_dimreg_pal[ match(cd4_pal, names_gdm_dimreg_pal )]


# Functions to create plots for panicos project

show_cluster <- function(x, cv, tt="", reduction="tsne", col="red"){
  ccells = WhichCells(x, expression=seurat_clusters==cv)
  DimPlot(x,cols.highlight = col, cells.highlight=ccells, reduction=reduction, sizes.highlight=0.5) + labs(title=tt) +
    NoLegend()
}

show_feature <- function(x, feature, reduction="tsne"){
  FeaturePlot(x, feature, reduction=reduction)
}

iter_show_feature <- function(x, features, reduction="tsne"){
  for(i in features){
    print(show_feature(x, i, reduction));
    invisible(readline(prompt=""))
  }
}

volcano_de_plot <- function(x, titles="", nsig=10){
  # Set significance
  x$sig = ifelse(x$p_val_adj <=0.05 & x$avg_log2FC >= 1, "DE", "Unsig")
  
  # Get top genes
  topg = dplyr::filter(x, sig=="DE") %>% dplyr::top_n(nsig, avg_log2FC) %>% dplyr::pull(gene)
  
  ggplot(x, aes(x=pct.1-pct.2, y=avg_log2FC, size=-log2(p_val_adj), color=sig)) + 
    geom_point(alpha=0.5) + 
    scale_x_continuous(limits=c(NA,1)) + 
    theme_classic() +
    geom_text_repel(aes(label=ifelse(sig=="DE", gene, "")), size=4, max.overlaps=50) +
    scale_color_manual(values=c("DE"="red", "Unsig"="grey")) +
    labs(title=titles)

}

#' Best volcano so far
volcano_de2 <- function(x){
  df_sig = dplyr::filter(x, abs(avg_log2FC)>=0.5 & p_val_adj <= 0.01)
  ggplot(x, aes(x=avg_log2FC, y=-log10(p_val_adj))) + geom_point(alpha=.5, size=1.2) + theme_bw() + 
    labs(x="avg_log2FC") + geom_hline(yintercept=-log10(0.01), linetype="dashed", color="grey") + geom_vline(xintercept = c(-.5,.5), linetype="dashed", color="grey") +
    geom_point(data=df_sig, color=ifelse(sign(df_sig$avg_log2FC)==-1, "blue", "red")) + 
    ggrepel::geom_text_repel(aes(label=gene), data=df_sig)
}

volcano_de_clust2 <- function(x, top_de){
  dc = ggplot(x, aes(x=avg_log2FC, y=-log10(p_val_adj))) + geom_point(alpha=.6, size=1, color="grey") + theme_bw() + 
    labs(x="avg_log2FC") + geom_hline(yintercept=-log10(0.01), linetype="dashed", color="grey") + geom_vline(xintercept = c(-.5,.5), linetype="dashed", color="grey")
  volc = dc + geom_point(data=top_de, aes(fill=factor(cluster)), color="black", size=3, shape=21) + ggrepel::geom_text_repel(aes(label=gene), data=top_de)
  return(volc)
  #geom_point(data=df_sig, color=ifelse(sign(df_sig$avg_log2FC)==-1, "blue", "red"))# + 
    #ggrepel::geom_text_repel(aes(label=gene), data=df_sig)
}

volcano_de_clust <- function(x){
  # Set top 10 DE genes by seurat cluster
  cluster_ids = as.character(unique(x$ident))
  
  cplots <- lapply(cluster_ids, function(i){
    volcano_de_plot(dplyr::filter(x, ident==i), i)
  })
  
  patchwork::wrap_plots(cplots)
  
}

cluster_by_case <- function(x, cluster){
  TSNEPlot(subset(x, seurat_clusters==cluster), split.by="control")
}

show_count_pat_clust <- function(x){
  ggplot(x@meta.data, aes(x=forcats::fct_reorder(multi_q, control=="CASE"), fill=control)) + 
    geom_bar(color="black", size=0.5) + scale_fill_viridis_d() + 
    facet_wrap(~seurat_clusters) + theme_classic() + 
    theme(axis.text.x=element_blank()) + labs(x="", y="N cells")
}

show_feat <- function(x, feat){
  FeaturePlot(x, feat, reduction="tsne") / VlnPlot(x, feat)
}

plot_med_per_pat_clust <- function(x){
  ggplot(x, aes(x=seurat_clusters, y=med)) + geom_point() + geom_errorbar(aes(ymin=med-(iqr/2), ymax=med+(iqr/2), width=0.15)) + 
    theme_bw() + labs(title=x$label[[1]], y="Medain number of cells per patient") + scale_y_continuous(n.breaks=20)
}


plot_cells_per_pat_median <- function(treg, cd4){
  plot_med_per_pat_clust(cells_per_pat_clust(treg, "treg")) /
    plot_med_per_pat_clust( cells_per_pat_clust(cd4, "cd4"))
}


show_dimplots <- function(x,y){
  px = (UMAPPlot(x, cols=gdm_dimreg_pal) + NoLegend()) | (TSNEPlot(x, cols=gdm_dimreg_pal))
  py = (UMAPPlot(y, cols=gdm_dimreg_pal) + NoLegend()) | (TSNEPlot(y, cols=gdm_dimreg_pal))
  px / py
}

show_tsneplot <- function(x, pal=gdm_dimreg_pal){
  plt = (TSNEPlot(x,  cols=pal)) + NoLegend() # + theme(legend.position="bottom")
  LabelClusters(plt, id="ident", box=T, color="white")
}


plot_cond_tsne <- function(x, col_gdm){
  Idents(x) <- "cond_cln"
  TSNEPlot(x, cols=c(col_gdm, gdm_dimreg_pal[12]))  + theme(legend.position="bottom") + labs(title="")
}

theme_cust_border <- function(){
  theme(panel.border=element_rect(fill=NULL, color="black"),
        axis.line=element_blank())
}

ShowSCTClusterDeTop <- function(x, de, reduction="tsne"){
  DefaultAssay(x) <- "SCT"
  cc = sort(unique(x$seurat_clusters))
  
  # Add color based on annotation
  de$color = ifelse(de$CONTROL_avg_log2FC >= 0, "red", "blue")
  
  # Rank de top 3 genes based on fold-change. All significant at this point.
  de3 <- dplyr::group_by(de, ident) %>% dplyr::top_n(5, abs(CONTROL_avg_log2FC))
  
  # Split de by cluster
  desplit <- split(de3, de3$ident)
  
  # Loop over cluster and generate plot for all three top genes
  deplots <- lapply(
    names(desplit),
    function(cluster){
      data = desplit[[cluster]]
      p_highlight = show_cluster(x, cluster, cluster, reduction=reduction)
      p_tops <- lapply(
        seq(nrow(data)),
        function(i){
          gene = data[i,]$gene
          col = data[i,]$color
          FeaturePlot(x, gene, reduction=reduction, cols=c("lightgrey", col))
        }
      )
      patchwork::wrap_plots(c(list(p_highlight), p_tops), ncol=4)
    })
  
  # Use patchwork to combine plots into single column plot
  p_topw <- patchwork::wrap_plots(deplots, ncol=1)
  return(p_topw)
}


# Plot volcano of clusters
volcano_within_clusters <- function(x, label, nrow=5, cutoff=0.25){
  
  p_volcano <- function(x){
    # Generate label dataset
    df_sig = dplyr::filter(x, abs(avg_log2FC)>=cutoff & p_val_adj <= 0.05)
    df_sig$color = ifelse(sign(df_sig$avg_log2FC)==-1, "down", "up")
    
    ggplot(x, aes(x=avg_log2FC, y=-log10(p_val_adj))) + 
      # Set style and show cut-offs
      geom_hline(yintercept=-log10(0.01), linetype="dashed", color="grey") +
      geom_vline(xintercept = c(-cutoff,cutoff), linetype="dashed", color="grey") +
      theme_bw() +
      # Plot all genes
      geom_point(alpha=.5) +
      # Label significant genes
      geom_point(data=df_sig, aes(color=color)) + 
      ggrepel::geom_text_repel(aes(label=gene), data=df_sig) +
      labs(x="avg_log2FC", title=paste0(unique(x$seurat_clusters), "_", unique(x$cluster))) +
      guides(color="none") +
      scale_color_manual(values=c("down"="blue", "up"="red"))
  }
  
  # Iterate over clusters and generate plots
  clust_iter <- sort(unique(x$seurat_clusters))
  plots <- mclapply(
    clust_iter, function(i){
      p_volcano(dplyr::filter(x, seurat_clusters==i))
    }
  )
  plt <- patchwork::wrap_plots(plots, nrow=nrow) + patchwork::plot_annotation(title=label)
  return(plt)
}


de_dot_plot <- function(x, col.scale=c("blue", "grey", "red")){
  # x = DE results table
  # Filter for significant genes
  x$is_sig = ifelse(x$is_sig==T & !is.na(x$is_sig), T, F)
  ggplot(x, aes(y=gene, x=cluster, fill=avg_log2FC, size=-log10(p_val_adj))) + 
    geom_point(aes(shape=is_sig, alpha=is_sig)) +
    #geom_point(aes(fill=avg_log2FC, size=-log10(p_val_adj)*0.75), shape=21, alpha=0.5) +
    #geom_point(aes(fill=avg_log2FC, size=-log10(p_val_adj)), shape=22, data=dplyr::filter(x, is_sig==T)) +
    theme_bw() +
    labs(x="", y="") +
    scale_x_discrete(position="bottom") +
    theme(axis.text.x=element_text(angle=45,hjust=1.05, vjust=1.05)) +
    scale_fill_gradient2(
      low=col.scale[[1]], mid=col.scale[[2]], high=col.scale[[3]],
      n.breaks=8) +
    scale_color_manual(
      values=c("FALSE"="white", "TRUE"="black"),
      breaks = c(TRUE, FALSE), # Values to include in the legend
      labels = c("q-value ≤ 0.05", "non-significant") # Labels to replace TRUE and FALSE
    ) +
    scale_alpha_discrete(range=c(0.55, 1), labels = c('N.S.', 'adjusted p ≤ 0.05')) +
    scale_shape_manual(values=c(21, 22), labels=c('N.S.', 'adjusted p ≤ 0.05')) +
    labs(shape="Significance", size="-log10\n(adjusted p)", fill="Average Log2\nFold-change", alpha='Significance')
}

subassign_plot <- function(x, features){
  (FeaturePlot(x, features, col=c("white", "blue"), reduction="tsne") | TSNEPlot(x, cols=gdm_dimreg_pal)) /
    VlnPlot(x, features, pt.size=0, cols=gdm_dimreg_pal)
}

subassign_thresh <- function(x, feature, thresh, features){
  x = subset(x, !!sym(feature)>=thresh)
  subassign_plot(x, features) + plot_annotation(title=paste0(feature, ">=", thresh))
}

subtype_heatmap <- function(avg_per_clust, genelist, is.treg=T){
  hm = DoHeatmap(avg_per_clust, size=5, group.bar=F,#group.by="clust_group",
                 features=genelist, angle=0, hjust=0.3,
                 assay="SCT", draw.lines=F,
                 group.colors = gdm_dimreg_pal[c(1,6,10)]) + guides(color="none") + 
                theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  
  if(is.treg){
    hm = hm +
      scale_fill_gradientn(colors=c(gdm_pal[["teal"]], gdm_pal[["white"]], gdm_pal[["orange"]]))
  } else {
    hm = hm +
      # Make consistent as colours for CD4 and Treg are not informative
      #scale_fill_gradientn(colors=c(gdm_pal[["darkgrey"]], gdm_pal[["white"]], gdm_pal[["pink"]]))  +
      scale_fill_gradientn(colors=c(gdm_pal[["teal"]], gdm_pal[["white"]], gdm_pal[["orange"]]))
  }
  return(hm)
}

show_condition_tsne <- function(x){
  TSNEPlot(x, group.by="cond_cln", cols=cond_cln_pal)
}


# Abundance Functions
plot_abundance_fbc <- function(x, tt, pal=gdm_dimreg_pal){
  ggplot(x@meta.data, aes(x=cond_cln)) + 
    geom_bar(aes(fill=cluster), position="fill", color="black") +
    scale_fill_manual(values=pal) + theme_classic() + 
    labs(x="", y="Proportion", fill="", title=tt) +
    theme(legend.position="bottom") + 
    theme(axis.text.x=element_text(angle=45, vjust=0.5))
}

plot_propeller_res <- function(x, tt, sig_col){
  ggplot(x, aes(y=cluster, x=PropRatio)) + 
    geom_vline(xintercept=c(.75, 1.25), linetype="dashed", color="grey") +
    geom_vline(xintercept=1, color="black") +
    geom_point(size=3.5) +
    geom_point(color=sig_col, size=2, data=x %>% dplyr::filter(abs(FDR)<=0.1)) +
    theme_bw() +
    guides(color="none") +
    labs(x="Proportion Ratio (GDM:Control)", y="", title="") +
    scale_x_continuous(n.breaks=6)
}

# Plot3: proportions
plot_pat_proportion <- function(x, col){
  case_levels <- x@meta.data %>% dplyr::select(multi_q, cond_cln) %>%
    dplyr::arrange(cond_cln, multi_q) %>% unique() %>% dplyr::pull(multi_q)
  
  pdata <- x@meta.data %>% 
    dplyr::group_by(multi_q) %>%
    dplyr::mutate(pat_cells = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(multi_q, cond_cln, cluster, pat_cells) %>%
    dplyr::summarise(clust_cells=n()) %>%
    dplyr::mutate(clust_prop=clust_cells/pat_cells) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(multi_q = fct_relevel(multi_q, case_levels))
  
  ggplot(pdata,aes(x=multi_q, y=clust_prop, fill=cond_cln)) + 
    geom_col() + theme_bw() +
    scale_y_continuous(n.breaks=3) +
    facet_wrap(~fct_rev(cluster), ncol=1, scales="free_y", strip.position="right") +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
    scale_fill_manual(values=c("black", col)) + 
    theme(strip.background = element_rect(fill="white", color=NULL),
          legend.position="top", # strip.text=element_text(size=4),
          axis.text.x=element_blank()) +
    labs(x="", fill="", y="Proportion") 
}

# Relevel factor to split by gdm patients
de_viol_plots <- function(x, obj, label, cols, figdir){
  DefaultAssay(obj) <- "SCT"
  pviol_files = mclapply(
    unique(x$seurat_clusters),
    function(i){
      print(sprintf("DE VIOL::%s", i))
      dsub = dplyr::filter(x, seurat_clusters == i & is_sig)  
      plt = VlnPlot(obj, dsub$gene, stack=T, flip=T, group.by="multi_q", split.by="cond_cln",
                    cols=cols) + labs(title=paste0(label, "_", dsub$cluster[[1]]))
      outfile = fs::path(figdir, paste0("vln_", label, "_", dsub$seurat_clusters[[1]], ".png"))
      ggsave(outfile, plt, width=8, height=round(
        length(dsub$gene)*1, 0))
    },
    mc.cores=6
  )
  return(pviol_files)
}

#gsea plots
plot_nes <- function(x, title){

ggplot(x, aes(
  x=forcats::fct_reorder(pathway, NES), y=NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=title) + 
  theme_minimal()

}

plot_enrich <- function(wed, gmts, pathway, add_title=T){
  gmtl <- Reduce(c, (lapply(gmts, gmtPathways)))
  stopifnot(pathway %in% names(gmtl))
  pathway_genes = gmtl[[pathway]]

  gene_ranks = wed$avg_log2FC
  names(gene_ranks) = wed$gene
  gene_ranks = gene_ranks[order(gene_ranks)]
  
  penr = plotEnrichment(pathway_genes, gene_ranks)
  
  if(add_title){
    penr = penr + labs(title=paste0(
    wed$dataset[[1]], "_", wed$cluster[[1]], " ::: ", pathway
  ))
  }
  
  return(penr)
}

plot_enrich_set <- function(gsea_res, de_res, gmts){
  plts <- lapply(
    seq(nrow(gsea_res)),
    function(i){
      # Get the current pathway
      cur_data = gsea_res[i,]
      pathway = cur_data$pathway
      # Filter de_res to the dataset and cluster
      cur_de = de_res %>% dplyr::filter(
        dataset==cur_data$dataset & seurat_clusters == cur_data$seurat_clusters
      )
      # Create plot
      plot_enrich(cur_de, gmts, pathway) + labs(
        title=paste0(cur_data$dataset, "_", cur_data$cluster, " ::: ", cur_data$pathway)
      )
      
    }
  )
  patchwork::wrap_plots(plts, ncol=1)
  
}


plot_nes_wrap <- function(clusters, res_gsea, prefix, color){
  lapply(
    unique(clusters),
    function(o){
      e = filter(res_gsea, cluster==o)
      plot_nes(e, prefix) + 
        scale_fill_manual(values=c("FALSE"="lightgray","TRUE"=color)) 
      
    }
  ) %>% patchwork::wrap_plots()
}

plot_found_wrap <- function(clusters, res_gsea, prefix, color, pathway_filter){
  lapply(
    unique(clusters),
    function(o){
      e = filter(res_gsea, cluster==o) %>% group_by(gmt) %>% top_n(10, abs(NES))
      plot_nes(e, paste0(prefix, o)) + 
        scale_y_continuous(limits=c(-3,3)) +
        scale_fill_manual(values=c("FALSE"="lightgray","TRUE"=color)) +
        facet_wrap(~gmt, nrow=1, scales="free_y")
      
    }
  ) %>% patchwork::wrap_plots(ncol=1)
}



build_escape_heatmap <- function(x, genesets, heatmap.colors){

  # Deprecated: inclust removed to be more flexible
  # # Downsample to N largest by selected clusters (excluding Other so not all cells shown)
  # count_selclust = table(x$inclust)
  # count_selclust = count_selclust[names(count_selclust)!="Other"]
  # nmax = count_selclust[which.max(count_selclust)]
  # Idents(x) = "inclust"
  # tsd = subset(x, downsample=nmax )
  # 
  # Filter any cells that are NA across any of the features
  tsd = filter_na_cells(x, genesets)
  
  # Plot heatmap
  colorblind_vector <- colorRampPalette(heatmap.colors)
  #png(filename=outfile, width=15, height=5, units="in", res=300)
  plt <- 
    ggplotify::as.ggplot(
      dittoHeatmap(tsd, genes = NULL, 
               metas = genesets, 
               heatmap.colors = colorblind_vector(50),
               annot.by = c("cond_cln"),
               annot.colors = rev(gdm_pal[c("white", "purple","green", "orange", "lightgrey", "pink", "teal", "ogreen", "orange")]),
               cluster_cols = T,
               cluster_rows=F,
               fontsize = 7)
    )
  return(plt)
  #dev.off()
}

p_wgcna_dend <- function(x, title, outfile){
  png(outfile, width=10, height=5, units="in", res=300)
  PlotDendrogram(conn, main=title)
  dev.off()
}

# Slingshot

plot_sling_lineage <- function(sling, lin, labels=T){
  # ss is a slingshot object
  # Default plots curves and tree, but may need to choose one for final pub.
  curves <- slingCurves(sling$ss, as.df=T)
  nlin = max(curves$Lineage)
  stopifnot(lin <= nlin)
  mst <- slingMST(sling$ss, as.df=T)
  
  fill_labs = ifelse(labels, "scc", "sc")
  gp <- ggplot(sling$ts, aes(x=tSNE_1, y=tSNE_2)) + 
    geom_point(aes(fill=!!sym(fill_labs)), alpha=.3, shape=21, color="transparent") + theme_bw() + 
    scale_fill_manual(values=gdm_dimreg_pal)
  
  gp + geom_point(data=mst, size=3) + 
    geom_path(data = curves  %>% dplyr::arrange(Order) %>% dplyr::filter(Lineage==lin), color="darkred") +
    geom_path(data=mst %>% dplyr::arrange(Order) %>% dplyr::filter(Lineage==lin), size=1, alpha=0.7, color="red") +
    labs(title=glue("Slingshot:{lin}"))
}


plot_lineage <- function(seu){
  lins = purrr::keep(names(seu@meta.data), ~stringr::str_detect(.x, "^Lineage"))
  if(length(lins)==0){
    log_info("No lineage present. Have you annotated pseudotime?")
    return(NA)
  }
  linplots = c(list(show_tsneplot(seu)), lapply(lins, function(e) FeaturePlot(seu, e, reduction="tsne")))
  wlinplots = patchwork::wrap_plots(linplots, ncol=2, guides = "collect") & theme(legend.position="bottom")
  return(wlinplots)
}

plot_lineage_vln <- function(seu){
  Idents(seu) = "cluster"
  lins = purrr::keep(names(seu@meta.data), ~stringr::str_detect(.x, "^Lineage"))
  if(length(lins)==0){
    log_info("No lineage present. Have you annotated pseudotime?")
    return(NA)
  }
  lin_viols = lapply(
    lins, function(n) VlnPlot(seu, n, col=gdm_dimreg_pal, pt.size = 0)
  )
  linwrap = wrap_plots(lin_viols, ncol=2, guides="collect") & theme(legend.position="bottom")
  return(linwrap)
}

plot_lin_cond_box <- function(seu){
  lins = purrr::keep(names(seu@meta.data), ~stringr::str_detect(.x, "^Lineage"))
  lclong = seu@meta.data %>% dplyr::select(
    seurat_clusters, cluster, cond_cln, starts_with("Lineage")
  ) %>% tidyr::pivot_longer(cols=starts_with("Lineage"), names_to="lin", values_to="pseu")
  ggplot(lclong, aes(x=lin, y=pseu, fill=cond_cln)) + geom_boxplot(outlier.shape=NA) + theme_bw() +
    scale_fill_manual(values=as.character(gdm_pal[c("grey", "orange")]))
}


plot_lin_cond_hist <- function(seu){
  lins = purrr::keep(names(seu@meta.data), ~stringr::str_detect(.x, "^Lineage"))
  lclong = seu@meta.data %>% dplyr::select(
    seurat_clusters, cluster, cond_cln, starts_with("Lineage")
  ) %>% tidyr::pivot_longer(cols=starts_with("Lineage"), names_to="lin", values_to="pseu")
  ggplot(lclong, aes(x=pseu, fill=cond_cln)) + geom_density(alpha=0.5) + theme_bw() +
    scale_fill_manual(values=as.character(gdm_pal[c("grey", "orange")])) + facet_wrap(~lin, ncol=1) +
    theme(strip.background = element_rect(fill="white"))
}

plot_gdm_meta <- function(data){
  p_ethnicity <- ggplot(data, aes(x = ethnicity, fill = status_c)) +
    geom_bar(position = "dodge", color="black", alpha=0.5) +
    labs(title = "", x = "Ethnicity", y = "Count") +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    theme_bw() +
    scale_y_continuous(n.breaks=10) +
    theme(legend.position="none") + labs(fill="", y="")
  
  p_bmi <- ggplot(data, aes(x = status_c, y=bmi, fill = status_c)) +
    geom_boxplot(alpha=0.5, outlier.shape=NA) +
    geom_jitter(alpha=0.5, width=0.13, shape=21) +
    labs(title = "", x = "", y = "Body mass index (BMI)") +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    scale_y_continuous(limits=c(min(data$bmi)-5,NA), n.breaks=8) +
    theme_bw() + theme(legend.position="none")

  p_map <-  ggplot(data, aes(x = status_c, y=map, fill = status_c)) +
    geom_boxplot(alpha=0.5, outlier.shape=NA) +
    geom_jitter(alpha=0.5, width=0.13, shape=21) +
    labs(title = "", x = "", y = "Mean arterial pressure (MAP)") +
    scale_fill_manual(values = c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    scale_y_continuous(limits=c(min(data$map)-5,NA), n.breaks=8) +
    theme_bw() + theme(legend.position="none")
  
  p_meta <- patchwork::wrap_plots(p_ethnicity, p_bmi, p_map, ncol=1)
  return(list(plot=p_meta, indiv=list(p_bmi, p_ethnicity, p_map)))
}

# Covariates.R plots

plot_stirms_box <- function(strims, figdir){
  bulk_plot = ggplot(strims, aes(x=status, y=expr)) + geom_boxplot(outlier.shape=NA) + facet_wrap(~gene)
  ggsave(fs::path(figdir, "strims_boxplots.png"), bulk_plot, width=10, height=10)
  
}

# Performs t-test and add significance stars
plot_stirms_gene <- function(strims){
  ggplot(strims, aes(x=status, y=expr, fill=status)) + geom_violin(alpha=0.8) +
    stat_summary(fun = "median",geom='point', shape='-', size=40, alpha=0.5) +
    geom_signif(comparisons=list(c("GDM", "NGT"))) + theme_bw() + theme(legend.position="none")+
    labs(title=strims$gene[1]) +
    scale_fill_manual(values=c(gdm_pal[["orange"]], gdm_pal[["lightgrey"]]))
}

# Plot for all genes in strims
p_wrap_stirms_gene <- function(strims, figdir){ 
  strims_genes = unique(strims$gene)
  strims_gene_plots <- lapply(
    unique(strims$gene),
    function(x){ plot_stirms_gene ( dplyr::filter(strims, gene == x) ) }
  )
  p_sgp <- patchwork::wrap_plots(strims_gene_plots, ncol=6)
  ggsave(fs::path(figdir, "strims_violins_pval.png"), p_sgp, width=20, height=20)
  return(p_sgp)
}

# CD4 covariates top
plot_cd4_top_genes_box <- function(dbox){
  dbox %>%
    ggplot(aes(x=gene, y=expr, fill=status_c)) + 
    geom_violin() +
    stat_summary(fun = median, fun.min=median, fun.max=median,
                 color=gdm_pal[["yellow"]], geom = "crossbar", width=0.3, position= position_dodge(width = 0.9)) +
    theme_bw() + labs(x="", y="Log10(CPM)") +
    scale_fill_manual(values=c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]])) +
    geom_jitter(shape=21,
                position=position_jitterdodge(dodge.width=0.9, jitter.width=0.05),
                alpha=0.5) +
    facet_wrap(~dataset, ncol=1) +
    theme_light()
}

# THrowaway plot for cd4 dimreg
plt_umap_dimreg <- function(cd4_dataset, gtop, neighbors=5, ptitle=""){
  pca_result <- prcomp(cd4_dataset[, gtop], scale = F)
  
  # Extract PCA scores
  pca_scores <- as.data.frame(pca_result$x[, 1:2])
  pca_scores$status <- ifelse(cd4_dataset$status==1, "GDM", "CONTROL")
  
  # Plot PC1 and PC2 using ggplot
  plt_pca = ggplot(pca_scores, aes(x = PC1, y = PC2, color = status)) +
    geom_point() +
    labs(title=ptitle, x = "Principal Component 1", y = "Principal Component 2") +
    theme_minimal() +
    scale_color_manual(values=c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]]))
  
  umap_result <- umap(cd4_dataset[, gtop], n_neighbors = neighbors)
  
  # Combine UMAP results with status column
  umap_data <- data.frame(cbind(umap_result$layout, data.frame(status = ifelse(cd4_dataset$status==1, "GDM", "CONTROL"))))
  
  # Create a scatter plot using ggplot
  plt_umap = ggplot(umap_data, aes(x = X1, y = X2, fill = status)) +
    geom_point(size=7, shape=21, color="black") +
    labs(title=ptitle, x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
    theme_minimal() +
    scale_fill_manual(values=c(gdm_pal[["lightgrey"]], gdm_pal[["orange"]]))
  
  plt_dimreg = patchwork::wrap_plots(plt_umap, plt_pca, nrow=1) + plot_layout(guides="collect")
  plt_dimreg
  return(list(plt_pca, plt_umap))
}

make_meta_coef <- function(meta_mod){
  # Plots the metadata coefficient barplot
  meta_scores = meta_mod %>% select(model_name, contains("all_y")) %>% unique() %>% arrange(-all_y_auc)
  ## Fit to all data and get coefficients
  comod = meta_mod %>% filter(model_name=="status;ethnicity;bmi;map")
  extract_info <- function(row, names_vec=c("intercept", "ethnicity_white", "bmi", "map")){
    csum =  coef(summary(row$model[[1]]))
    csum = data.frame(t(csum))
    names(csum) = names_vec
    csum$feats = c("estimate", "std_err", "z_value", "p_value")
    cf = cbind(row %>% select(model_name, fold_no), csum)
    return(cf)
  }
  meta_fold = rbindlist(lapply(seq(nrow(comod)), function(i) extract_info(comod[i,])))
  meta_foldl =  meta_fold  %>% tidyr::pivot_longer(cols=c(-model_name, -fold_no, -feats)) %>% dplyr::mutate(is_intercept=ifelse(name=="intercept", "Intercept", "Features")) %>% mutate(dataset="combined")
  metae= dplyr::filter(meta_foldl, feats=="estimate")
  metap = dplyr::filter(meta_foldl, feats=="p_value")
  
  
  comod2 = meta_mod %>% filter(model_name=="status;bmi")
  meta_fold2 = rbindlist(lapply(seq(nrow(comod2)), function(i) extract_info(comod2[i,], c("intercept", "bmi"))))
  meta_foldl2 =  meta_fold2  %>% tidyr::pivot_longer(cols=c(-model_name, -fold_no, -feats)) %>% dplyr::mutate(is_intercept=ifelse(name=="intercept", "Intercept", "Features")) %>% mutate(dataset="bmi_only")
  metae2 = dplyr::filter(meta_foldl2, feats=="estimate")
  metap2 = dplyr::filter(meta_foldl2, feats=="p_value")
  
  meta_plot_e = rbind(meta_foldl, meta_foldl2) %>% dplyr::filter(feats %in% c("estimate", "p_value")) %>% 
    tidyr::pivot_wider(names_from='feats', values_from='value')
  
  p_meta_coef = ggplot(meta_plot_e, aes(x=name, y=estimate, group=as.factor(fold_no))) + 
    geom_point(aes(
      y=(0.75*sign(estimate)),
      color=ifelse(p_value<=0.1, "blue", "white")
    ),
    position=position_dodge(width=0.9)) +
    geom_point(aes(
      y=(0.75*sign(estimate)), 
      color=ifelse(p_value<=0.05, "grey", "white")
    ),position=position_dodge(width=0.9), shape=5) +
    geom_col(position="dodge", color='black', fill=gdm_pal[['orange']]) + theme_bw() +
    labs(x="Feature", y="Coefficient") + ggh4x::facet_grid2(dataset~is_intercept, scales = "free", independent="y") +
    theme(legend.position="none") + scale_color_manual(values=c(gdm_pal[["orange"]], "red", "white"))
  return(p_meta_coef)
}

plot_roc_data <- function(roc_data){
  
  # Plot ROC curves
  roc_text = dplyr::select(roc_data, dataset, model_name_strims, all_y_auc) %>% unique() %>%
    dplyr::mutate(label=paste0(dataset, "-AUC:", round(all_y_auc, 2))) %>%
    dplyr::mutate(x=0.7, y=ifelse(dataset=="Strim", 0.2, 0.1))
  p_strim_roc = ggplot(roc_data %>% filter(dataset!="CD4+_Strim"), aes(m = y_prob, d = y_truth, color=dataset)) + 
    geom_roc(n.cuts=0, alpha=0.5) +
    geom_text(data=roc_text, mapping=(aes(x=x, y=y, label=label, color=dataset)), inherit.aes=F) +
    facet_wrap(~model_name_strims, ncol=5) + theme_bw() +
    theme(legend.position = "none") +
    scale_colour_manual(values=c(gdm_pal[["orange"]], gdm_pal[["lightgrey"]])) +
    labs(x="False Positive Fraction", y="True Positive Fraction")
  return(p_strim_roc)
}

plot_roc_data_single <- function(roc_data_ss, rgenes, ncol=2){
  rss = dplyr::filter(roc_data_ss, model_name %in% rgenes)
  udatasets = unique(rss$dataset)
  # Plot ROC curves
  roc_text = dplyr::select(rss, dataset, model_name, all_y_auc) %>% unique() %>%
    dplyr::mutate(label=paste0(dataset, "-AUC:", round(all_y_auc, 2))) %>%
    dplyr::mutate(y=dplyr::case_when(
      dataset==udatasets[4] ~ 0.4,
      dataset==udatasets[1] ~ 0.3,
      dataset==udatasets[2] ~ 0.2,
      dataset==udatasets[3] ~ 0.1), x=0.7)
  p_strim_roc = ggplot(rss, aes(m = y_prob, d = y_truth, color=dataset)) + 
    geom_roc(n.cuts=0, alpha=0.5) +
    geom_text(data=roc_text, mapping=(aes(x=x, y=y, label=label, color=dataset)), inherit.aes=F) +
    facet_wrap(~model_name, ncol=ncol) + theme_bw() +
    theme(legend.position = "none") +
    scale_colour_manual(values=c(gdm_pal[["lightgrey"]], gdm_pal[["teal"]], gdm_pal[["pink"]], gdm_pal[["orange"]])) +
    labs(x="False Positive Fraction", y="True Positive Fraction") +
    theme(strip.background = element_rect(fill="white"))
  return(p_strim_roc)
}

plot_roc_f1_single <- function(x){
  df = x %>% dplyr::select(-fold_no, -fold_id, -y_truth, -y_prob) %>% 
      unique() #%>% dplyr::filter(all_y_auc >= 0.5 & all_y_f1 >=0.5)
  df <- df %>%
    mutate(rank_metric = all_y_f1 + all_y_auc) %>%
    arrange(desc(rank_metric))
  df$is_pseu = ifelse(
    stringr::str_detect(df$dataname, "scRNA"),
    "Single-cell Pseudobulk",
    "Bulk RNA"
  )
  datasets = split(df, df$is_pseu)
  
  plt <- ggplot(datasets[[1]], aes(x = all_y_f1, y = all_y_auc)) +
    geom_line(aes(group=model_name), linetype='dashed', alpha=0.6) +
    geom_jitter(aes(color=dataname), size=3, alpha=.5, shape=17) +
    geom_label_repel(aes(label=model_name), max.overlaps=50,
                     alpha=0.8, box.padding = 0.25, size=3.5, force=10, max.iter=200) +
    theme_bw() +
    scale_color_manual(values=c(as.character(gdm_pal[c('green', 'orange')]))) +
    labs(x = "F1 Score", y = "AUC Value", color="Bulk RNA") +
    theme(legend.position=c(0.75,0.15), legend.box.background = element_rect(color="black", size=0.5)) +
    scale_x_continuous(limits=c(0.25,1)) +
    scale_y_continuous(limits=c(0.25,1)) 
  
  
  plt2 <- ggplot(datasets[[2]], aes(x = all_y_f1, y = all_y_auc)) +
    geom_line(aes(group=model_name), linetype='dashed', alpha=0.6) +
    geom_jitter(aes(color=dataname), size=3, alpha=.5) +
    geom_label_repel(aes(label=model_name), max.overlaps=50,
                     alpha=0.8, box.padding = 0.25, size=3.5) +
    theme_bw() +
    scale_color_manual(values=c(as.character(gdm_pal[c('darkgrey', 'pink')]))) +
    labs(x = "F1 Score", y = "AUC Value", color="Pseudobulk") +
    theme(legend.position=c(0.75,0.15), legend.box.background = element_rect(color="black", size=0.5)) +
    scale_x_continuous(limits=c(0.25,1)) +
    scale_y_continuous(limits=c(0.25,1)) 
  
  return(list(plt,plt2))
}

mtrnr_map_plot <- function(bp_pred){
ggplot(bp_pred, aes(x=map, y=expr, color=status_c)) + geom_point() +  scale_color_manual(values=c(gdm_pal[["orange"]], gdm_pal[["lightgrey"]])) + 
  geom_line(aes(y=pred)) +
  geom_text(data = dplyr::filter(bp_pred, expr>4.5), mapping=aes(label=id, x=map, y=expr), inherit.aes=F, vjust=-.2) +
  facet_wrap(~gene) +
  theme_light() +
  labs(x="Mean Arterial Pressure (MAP)", y="Expression [?]")
}

p_hallmark_founder_sig <- function(tssig, title){
  ggplot(tssig, aes(x = cluster, y = gene_set,
                    size=-log10(padj), fill=median.GDM-median.CONTROL,
                    color=ifelse(padj <=0.01, T, F))) +
    geom_point(shape=21) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="bottom") +
    scale_fill_gradient2(
      low=gdm_pal[["green"]],
      #mid = gdm_pal[["lightgrey"]],
      high = gdm_pal[["orange"]]
    ) +
    scale_color_manual(
      values=c("lightgrey", "black"),
      labels=c("Non-significant", "Adjusted p-value ≤ 0.01")
    ) +
    labs(x="", y="", title=title,
         color="", fill="Median GDM NES - Control NES")
}

plot_odds_rat <- function(x){
  # Create a ggplot for odds ratios and CIs
  ggplot(x, aes(x = odds_ratio, y = clean, xmin = lower_ci, xmax = upper_ci)) +
    geom_point(stat = "identity", position = position_dodge(width = 0.75), size=3) +
    geom_errorbar(position = position_dodge(width = 0.75), width = 0.2) +
    labs(
      x = "Odds Ratio",
      y = ""
    ) +
    geom_vline(xintercept=1, linetype="dashed", alpha=0.5)+
    theme_bw() +
    theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))
  
}

tops5_to_upset <- function(x){
  cup = x %>% dplyr::select(all_y_auc, model_name) 
  ugenes = unique(unlist((stringr::str_split(cup$model_name, ";"))))
  ugenes = ugenes[ugenes!="status"]
  ugenes_bool = sapply(ugenes, function(x) stringr::str_detect(cup$model_name, x))
  cup = cbind(cup, ugenes_bool)
  cup$status=NULL
  cupl = tidyr::pivot_longer(cup, cols=all_of(ugenes), names_to="gene", values_to="present")
  dcup = cupl %>% filter(present==TRUE) %>% group_by(model_name, all_y_auc) %>% summarise(genes=list(gene)) %>% dplyr::arrange(-all_y_auc)
  return(dcup)
}
