# Functions to run DE analysis

library(parallel)
library(Seurat)
library(ggrepel)
library(dplyr)
library(fgsea)
library(data.table)
library(MAST)
library(scProportionTest)

# Objects

# Extned sc_utils to enable plot storage
# Extend the perms objects
setOldClass(c("gg", "ggplot"))
setClass(
  "mySPT",
  contains=c("sc_utils", "ggplot"),
  slots=c(plot="ggplot")
) -> mySPT


cluster_de <- function(obj, label){
  DefaultAssay(obj) = "RNA"
  Idents(obj) = "seurat_clusters"
  cluster_ids = as.character(levels(obj))
  cluster_de = mclapply(cluster_ids, function(x){
    cl = FindMarkers(obj, ident.1 = x, only.pos=F)
    cl$gene = rownames(cl)
    cl$ident = x
    return(cl)
    },
    mc.cores=10
  )
  de_result = Reduce(rbind, cluster_de)
  de_result$label = label
  return(tibble(de_result))
}

de_top_filter <- function(x, type="all", fc.cutoff=0.5, top_n=2) {
  # Type can be pos, neg or all
  stopifnot(type %in% c("pos", "neg", "all"))
  
  if(type=="pos"){x = dplyr::filter(x, avg_log2FC>=fc.cutoff)}
  if(type=="neg"){x = dplyr::filter(x, avg_log2FC<=fc.cutoff)}
  
  x %>% dplyr::filter(
  abs(pct.2 - pct.1) >= 0.2 & # Minimum percentage
  p_val_adj <= 0.01) %>%
  dplyr::top_n(top_n, abs(avg_log2FC))
}

ShowClusterDeTop <- function(x, de, reduction="tsne"){
  DefaultAssay(x) <- "RNA"
  cc = sort(unique(x$seurat_clusters))
  
  get_cluster_top_genes <- function(x, de, i){
    pos_genes = de %>% dplyr::filter(ident==as.character(i)) %>% de_top_filter("pos")
    neg_genes = de %>% dplyr::filter(ident==as.character(i)) %>% de_top_filter("neg")
    fps = lapply(pos_genes$gene, function(y) FeaturePlot(x, y, reduction=reduction))
    if(length(fps)>2){fps = fps[1:2]}
    if(length(fps)==1){fps = list(fps[[1]], ggplot())}
    # Switch to under-expressed genes if none over-expressed
    if(length(fps)==0){
      fps = lapply(neg_genes$gene, function(y) FeaturePlot(x, y, reduction=reduction, cols=c("lightgray","red")))
      if(length(fps)>2){fps = fps[1:2]}
      if(length(fps)==1){fps = list(fps[[1]], ggplot())}
      if(length(fps)==0){fps = list(ggplot(), ggplot())}
    }
    plts = c(list(show_cluster(x, i, i, reduction=reduction)), fps)
    patchwork::wrap_plots(plts, ncol=3)
  }

  p_top <- lapply(cc, get_cluster_top_genes, x=x, de=de)
  p_topw <- patchwork::wrap_plots(p_top, ncol=1)
  return(p_topw)
}

call_de <- function(obj, idents, x.1, x.2){
  Idents(obj) = idents

  cl = FindMarkers(obj, ident.1=x.1, ident.2=x.2, logfc.threshold=0, test.use="MAST")
  cl$gene = rownames(cl)
  cl$ident = idents
  cl$x.1 = x.1
  cl$x.2 = x.2

  return(tibble(cl))
}

condition_de <- function(obj, label){
  Idents(obj) = "condition"
  conds = unique(obj$condition)
  stopifnot("CONTROL" %in% conds)
  conds = setdiff(conds, "CONTROL")
  
  cond_de = mclapply(conds, function(x){
    cl = FindMarkers(obj, ident.1=x, ident.2="CONTROL", logfc.threshold=0)
    cl$gene = rownames(cl)
    cl$ident = x
    cl
  }, mc.cores=10)
  
  de_result = Reduce(rbind, cond_de)
  de_result$label = label
  
  #de_result = dplyr::filter(de_result, p_val_adj <= 0.05)
  
  return(tibble(de_result))
}

condition_gsea <- function(
    de,
    gmt_f = fs::dir_ls(here::here("data/gsea_gene_sets/"), glob="*.gmt"),
    cores=6
    ){
  
  de_split = split(de, de$ident)
  
  gsea = mclapply(
    gmt_f,
    function(x){
      rbindlist(
        lapply(de_split, run_gsea, gmt=x),
        fill=T
      )
    },
    mc.cores=cores
  )
  
  gsea = rbindlist(gsea, fill=T)
  
  if(nrow(gsea)>0){
    # Filter if not null
    gsea = dplyr::filter(gsea, padj <= 0.1) # Soft filter
  }
  
  return(gsea)
}

run_gsea <- function(tbl, gmt, field="avg_log2FC"){
  # Ref for gmt_files :: fs::dir_ls(here::here("data/gsea_gene_sets/"), glob="*.gmt")

  # Run GSEA
  pathways = gmtPathways(gmt)
  ranks = tbl[[field]]
  names(ranks) = tbl$gene
  
  fgseaRes <- tryCatch(
    fgsea(pathways, ranks, minSize=3, maxSize=500),
    error=function(e) data.table()) # Min genes in set
  
  if(nrow(fgseaRes) > 0){
  fgseaRes$gmt = fs::path_file(gmt)
  fgseaRes = dplyr::arrange(fgseaRes, padj)
  #fgseaRes$itblnt = unique(i$itblnt)
  
  # Add metadata from input DE table
  fgseaRes$ident = tbl$ident[[1]]
  fgseaRes$label = tbl$label[[1]]
  }
  
  return(fgseaRes)
}

run_sc_prop_test <- function(x, label){
  # Create object for test:
  perm <- sc_utils(x)
  
  # Run permutation test 
  res = permutation_test(
    perm, cluster_identity = "cluster",
    sample_1 = "CONTROL", sample_2 = "CASE",
    sample_identity = "control"
  )
  res = as(res, "mySPT")
  res@plot = permutation_plot(res) + labs(title=paste0(label, "::"))
  
  # Run permutation test for each condition against the control
  # conditions <- unique(x@meta.data$control) %>% purrr::discard(~.x=="CONTROL")
  # perms <- mclapply(conditions, function(i){
  #   res = permutation_test(
  #     perm, cluster_identity = "cluster",
  #     sample_1 = "CONTROL", sample_2 = i,
  #     sample_identity = "condition"
  #   )
  #   res = as(res, "mySPT")
  #   res@plot = permutation_plot(res) + labs(title=paste0(label, "::", i))
  #   return(res)
  # }, mc.cores=4
  # )
  # names(perms) <- conditions

  return(res)
}

# Diffexp by condition and by cluster
run_de_cond_clust <- function(x, min_cells=30) {
  stopifnot(
    "condition" %in% names(x@meta.data) &
      "seurat_clusters" %in% names(x@meta.data)
  )
  
  Idents(x) = "seurat_clusters"
  # Filter any clusters rthat have too few cells
  combs = dplyr::count(x@meta.data, condition, seurat_clusters) %>% dplyr::filter(n>=min_cells)
  combs = split(combs, combs$condition == "CONTROL")
  # Select unique combo of condition and clusters to compare only if they are present in control
  indata = dplyr::filter(combs[["FALSE"]], seurat_clusters %in% combs[["TRUE"]]$seurat_clusters)
  
  # Define function for DE on a single condition and cluster
  # Here, we follow reocmmendation to use no logFC threshold and rank on this value
  ## https://github.com/ctlab/fgsea/issues/50#issuecomment-514599752
  de_cond_clust <- function(x, i, indata){
    cond = indata[[i, "control"]]; clust= indata[[i, "seurat_clusters"]]
    marks = FindMarkers(x, subset.ident=clust, group.by="control",
                        ident.1="CASE", ident.2="CONTROL", logfc.threshold = 0,
                        test.use="MAST")
    if(nrow(marks) > 0){
      marks$gene = rownames(marks)
      marks$condition = cond
      marks$cluster = clust
      marks$ident = paste(marks$condition, marks$cluster, sep="_")
      #marks = dplyr::filter(marks, p_val_adj <= 0.05)
    }
    return(marks)
  }
  
  res = rbindlist(
    mclapply(seq(nrow(indata)), function(i){
      de_cond_clust(x, i, indata)
    },
    mc.cores=5),
    fill=T)
  
  return(res)
}

call_de_cc <- function(x, min_cells=30) {
  stopifnot(
    "control" %in% names(x@meta.data) &
      "seurat_clusters" %in% names(x@meta.data)
  )
  
  DefaultAssay(x) = "RNA"
  Idents(x) = "seurat_clusters"
  # Filter any clusters that have too few cells
  combs = dplyr::count(x@meta.data, control, seurat_clusters) %>% dplyr::filter(n>=min_cells)
  combs = split(combs, combs$control == "CONTROL")
  # Select unique combo of condition and clusters to compare only if they are present in control
  indata = dplyr::filter(combs[["FALSE"]], seurat_clusters %in% combs[["TRUE"]]$seurat_clusters)
  
  
  # Define function for DE on a single condition and cluster
  # Here, we follow reocmmendation to use no logFC threshold and rank on this value
  ## https://github.com/ctlab/fgsea/issues/50#issuecomment-514599752
  de_cond_clust <- function(x, i, indata){
    cond = indata[[i, "control"]]; clust= indata[[i, "seurat_clusters"]]
    marks = FindMarkers(x, subset.ident=clust, group.by="control",
                        ident.1="CASE", ident.2="CONTROL", logfc.threshold = 0)
    if(nrow(marks) > 0){
      marks$gene = rownames(marks)
      marks$condition = cond
      marks$cluster = clust
      marks$ident = paste(marks$condition, marks$cluster, sep="_")
      #marks = dplyr::filter(marks, p_val_adj <= 0.05)
    }
    return(marks)
  }
  
  res = rbindlist(
    mclapply(seq(nrow(indata)), function(i){
      de_cond_clust(x, i, indata)
    },
    mc.cores=5),
    fill=T)
  
  return(res)
}

# Plots ----

# Plot volcano of clusters
volcano_within_clusters <- function(x, label, nrow=5){

  p_volcano <- function(x){
  # Generate label dataset
  df_sig = dplyr::filter(x, abs(avg_log2FC)>=0.5 & p_val_adj <= 0.01)
  df_sig$color = ifelse(sign(df_sig$avg_log2FC)==-1, "down", "up")
    
  ggplot(x, aes(x=avg_log2FC, y=-log10(p_val_adj))) + 
    # Set style and show cut-offs
    geom_hline(yintercept=-log10(0.01), linetype="dashed", color="grey") +
    geom_vline(xintercept = c(-.5,.5), linetype="dashed", color="grey") +
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


de_dot_plot <- function(x, col.scale){
  # x = DE results table
  # Filter for significant genes
  ggplot(x, aes(y=gene, x=cluster)) + 
    geom_point(aes(fill=avg_log2FC, size=-log10(p_val_adj)), shape=21) +
    theme_bw() +
    scale_x_discrete(
      position="top"
    ) + theme(
      axis.text.x = element_text(angle=45, vjust=-.1,  hjust=-.1)
    ) +
    labs(x="", y="") +
    scale_fill_gradient2(
      low=col.scale[[1]], mid=col.scale[[2]], high=col.scale[[3]],
      n.breaks=8)
}
