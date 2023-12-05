library(ggplot2)
library(dplyr)
library(forcats)
library(scProportionTest)
source("code/plots.R")

#' Run single cell proportion testing
#' @param x Seurat object
#' @param y label
sc_prop_test <- function(x, label, sample_identity, sample_1, sample_2, cluster_identity){
  # Create object for test:
  perm <- sc_utils(x)
  
  # Run permutation test 
  res = permutation_test(
    perm, cluster_identity = cluster_identity,
    sample_1 = sample_1, sample_2 = sample_2,
    sample_identity = sample_identity
  )
  return(res)
}

# Plot results from sc_prop_test
sc_prop_plot <- function(sco, xlabel="", title="", sig_col="red"){
  # sco is an object returned by running scProportionTest
  df = sco@results$permutation
  ggplot(df, aes(y=clusters, x=obs_log2FD)) + 
    geom_vline(xintercept=c(-0.5, 0.5), linetype="dashed", color="grey") +
    geom_vline(xintercept=0, color="black") +
    geom_pointrange(aes(xmin=boot_CI_2.5, xmax=boot_CI_97.5), color="black") +
    geom_point(color=sig_col, size=2, data=df %>% dplyr::filter(abs(obs_log2FD)>=0.5)) +
    theme_bw() +
    guides(color="none") +
    labs(x=xlabel, y="", title=title) +
    scale_x_continuous(n.breaks=5) +
    theme(axis.text.x=element_text(size=5))
}

# Plot cells per patient
sc_prop_pat <- function(obj, case.col="red"){
  # Error if expected fields are missing
  stopifnot(
    all(
      c("multi_q", "cond_cln", "cluster") %in% names(obj@meta.data)
    )
  )
  
  # Subset data
  pdata <- obj@meta.data %>%
    dplyr::select(multi_q, cond_cln, cluster)
  
  # Set row order for patient labels, ordering
  pdata <- pdata %>% 
    dplyr::mutate(
      multi_q = forcats::fct_reorder(
        multi_q,  ifelse(cond_cln=="CONTROL", 1, 2)
      )
    )
    
  # Plot
  ggplot(pdata,aes(x=multi_q, fill=cond_cln)) + 
    geom_bar() + theme_bw() +
      facet_wrap(~fct_rev(cluster), ncol=1, scales="free_y", strip.position="right") +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
    scale_fill_manual(values=c("black", case.col)) +
    theme(strip.background = element_rect(fill="white", color=NULL),
          legend.position="top", strip.text=element_text(size=4),
          axis.text.x=element_text(size=4)) +
    labs(x="", fill="")
  
}

# Calculate proportion of cells per condition field
sc_prop_cond <- function(x, group.by, color.by){
  meta = x@meta.data[, c(group.by, color.by)]
  ct = dplyr::count(meta, !!sym(group.by), !!sym(color.by))
  names(ct) = c("group", "cluster", "count")
  
  plt = ggplot(ct, aes(x=group, y=count)) + 
    geom_col(aes(fill=cluster), position="fill", color="black") +
    scale_fill_manual(values=gdm_dimreg_pal) + theme_classic() + 
    labs(x="", y="Cell Proportion", fill="") +
    theme(legend.position="bottom")
  
  # Return plot and counts
  return(list(plt, ct))
}
