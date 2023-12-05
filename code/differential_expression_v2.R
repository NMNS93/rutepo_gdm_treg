library(parallel)
library(Seurat)
library(MAST)
library(dplyr)
library(data.table)

call_de_cc <- function(x, min_cells=20) {
  stopifnot(
    "control" %in% names(x@meta.data) &
      "seurat_clusters" %in% names(x@meta.data)
  )
  
  DefaultAssay(x) = "SCT"
  Idents(x) = "seurat_clusters"
  

  #b.interferon.response <- FindMarkers(immune.combined.sct, assay = "SCT", ident.1 = "B_STIM", ident.2 = "B_CTRL",verbose = FALSE)
  # Filter any clusters that have too few cells
  combs = dplyr::count(x@meta.data, control, seurat_clusters) %>% dplyr::filter(n>=min_cells)
  combs = split(combs, combs$control == "CONTROL")
  # Select unique combo of condition and clusters to compare only if they are present in control
  indata = dplyr::filter(combs[["FALSE"]], seurat_clusters %in% combs[["TRUE"]]$seurat_clusters)
  
  # Define function for DE on a single condition and cluster
  # Here, we follow the recommendation to use no logFC threshold and rank on this value
  ## https://github.com/ctlab/fgsea/issues/50#issuecomment-514599752
  de_cond_clust <- function(x, i, indata){
    cond = indata[[i, "control"]]; clust=indata[[i, "seurat_clusters"]]
    # FindMarkers should have been performed previously for subtyping
    xs <- PrepSCTFindMarkers(subset(x, seurat_clusters==clust))
    marks = FindMarkers(xs, group.by="control",
                        ident.1="CASE", ident.2="CONTROL", logfc.threshold = 0, test.use="MAST")
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


call_de_all <- function(x, min_cells=20) {
  # Call differential expression between all clusters and gene sets
  # Follows: https://satijalab.org/seurat/articles/sctransform_v2_vignette
  
  # Check
  stopifnot(
    "control" %in% names(x@meta.data) &
      "seurat_clusters" %in% names(x@meta.data)
  )
  
  # Setup
  DefaultAssay(x) = "SCT"
  Idents(x) = "seurat_clusters"
  
  # Identify and filter any clusters with too few cells for DE in either case or control
  f_filter_clusters <- function(x, case_name){
    x@meta.data %>% dplyr::filter(control==case_name) %>% dplyr::count(seurat_clusters) %>% 
      dplyr::filter(n>=min_cells) %>% pull(seurat_clusters) %>% as.character %>% as.numeric()
  }
  pass_in_case = f_filter_clusters(x, 'CASE')
  pass_in_control = f_filter_clusters(x, 'CONTROL')
  analysis_clusters = intersect(pass_in_case, pass_in_control)
  
  # Define function for DE on a single cluster
  de_clust <- function(x, i){
    # x is seurat object
    # i is cluster identifier in seurat_clusters
    xsubset <- subset(x, seurat_clusters==i)
    
    xmarkers <- FindMarkers(
      xsubset, assay = "SCT", group.by='control', ident.1 = "CASE", ident.2 = "CONTROL",
      logfc.threshold=0, verbose=FALSE, recorrect_umi = FALSE, test.use='MAST', min.pct=0)

    if(nrow(xmarkers) > 0){
      xmarkers$gene = rownames(xmarkers)
      xmarkers$condition = 'CASE:CONTROL'
      xmarkers$seurat_clusters = factor(i)
      xmarkers$ident = paste(xmarkers$condition, xmarkers$cluster, sep="_")
      xmarkers <- dplyr::left_join(
        xmarkers,
        unique(x@meta.data[, c("label", "seurat_clusters", "clust_group", "cluster")])
      )
      xmarkers = data.table(xmarkers)
    }
    
    return(xmarkers)
  }
  
  
  # Apply DE for each cluster
  de_res = data.table::rbindlist(
    mclapply(
      analysis_clusters,
      function(e){ de_clust(x, e) },
      mc.cores=5
    )
  )
  
  de_res = data.table(de_res)
  
  return(de_res)
}

get_sig_de_genes_v3 <- function(x){
  # Identify significantly differentially expressed genes using the criteria:
  # - p_adj <= 0.05
  # - avg_log2FC (absolute) > 0.2
  # - pct expressed in both groups > 0.1
  
  # Applied after call DE pipeline to dataset
  x = data.table(x)
  signif = x[
    p_val_adj <= 0.05 & 
      (pct.1 >= 0.1 & pct.2 >=0.1) &
      abs(avg_log2FC) >= 0.2,
    
    .(gene, seurat_clusters)
  ]
  signif$is_signif = T
  
  return(signif)
  
}


get_sig_de_genes <- function(x){
#WARNING : Only returns top 10
  # Get top N differentially expressed genes.
  # Adds label `is_sig` if significant.
  
  # min.pct from DE FindMarkers function applies to either not both,
  # so we filter out spurious results using both
  x <- dplyr::filter(x, (pct.1>0.2 & pct.2>0.2))

  # Add a label for significant sites
  x <- dplyr::mutate(x, is_sig= abs(avg_log2FC)>.25 & p_val_adj <=.05)
  
  # Select top N by log2FC
  sig <- dplyr::group_by(x, seurat_clusters) %>%
      dplyr::top_n(10, abs(avg_log2FC))
  
  return(sig)
}

cluster_filter <- function(x){
  # Return list of clusters that have >3 patients with >20 cells in both patient groups
  dt <- x@meta.data[, c("cluster", "multi_q", "cond_cln")]
  # Get average cells per control patient per cluster to help set threshold for cells in min patients
  cpat_av = dt %>% dplyr::filter(cond_cln=="CONTROL") %>% dplyr::count(multi_q, cluster) %>% dplyr::summarise(mean(n)) %>% unlist()
  cthresh = cpat_av*0.25 # This is the threshold used for clusters
  cnts = dplyr::count(dt, cluster, multi_q, cond_cln) 
  clusts = dplyr::mutate(cnts, has_count = n > cthresh) %>%
    dplyr::group_by(cluster, cond_cln) %>%
    dplyr::summarise(n_count_pass = sum(has_count)) %>%
    dplyr::summarise(n_group_pass = all(n_count_pass>3)) %>%
    dplyr::filter(n_group_pass) %>%
    dplyr::pull(cluster)
  return(list(
    cpass=clusts,
    "threshold"=cthresh
  ))
}

order_patients_by_level <- function(obj){
  obj@meta.data$multi_q = as.character(obj@meta.data$multi_q)
  multi_q_order <- dplyr::select(obj@meta.data, multi_q, cond_cln) %>% dplyr::arrange(cond_cln) %>%
    dplyr::pull(multi_q) %>% unique()
  obj@meta.data$multi_q = forcats::fct_relevel(
    obj@meta.data$multi_q, multi_q_order
  )
  return(obj)
}

