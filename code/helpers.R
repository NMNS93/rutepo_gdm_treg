library(Seurat)

n_cells <- function(x) length(Cells(x))
med_umi_pc <- function(x) median(x@meta.data$nCount_RNA)
med_genes_pc <- function(x) median(x@meta.data$nFeature_RNA)
total_genes <- function(x) nrow(x)

n_gdm <- function(x) sum(x$control=="CASE")

mean_sd_cells <- function(x){
  ccounts = table(x@meta.data$multi_q)
  res = c(mean(ccounts, na.rm=T), sd(ccounts, na.rm=T))
  round(res, 1)
}

get_qc_stats <- function(x, label=NA){
  # Returns stats from object including:
  # number of cells, cells per condition, cells per patient
  tibble(
    cells = n_cells(x),
    avpat = avg_pat_cells(x),
    sdpat = avg_pat_cells(x, sd=T),
    cases = n_cases(x),
    controls = n_cases(x, case=F),
    cells_cases = n_cells_cases(x),
    cells_controls = n_cells_cases(x, case=F),
    label=label
  )
}

# QC stats helpers ----
n_cells <- function(x) length(Cells(x))
n_cases <- function(x, case=T){
  # Return n cells per case or control
  case_count = dplyr::select(x@meta.data, multi_q, control) %>% unique() %>% dplyr::count(control)
  if(case){
    return(case_count[case_count$control=="CASE",]$n)
  }else{
    return(case_count[case_count$control=="CONTROL",]$n)
  }
}
n_cells_cases <- function(x, case=T){
  if(case){
    n_cells(subset(x, control=="CASE"))
  }else{
    n_cells(subset(x, control=="CONTROL"))
  }
}
avg_pat_cells = function(x, sd=F){
  avpc = dplyr::count(x@meta.data, multi_q) %>% summarise(mm=mean(n), sd=sd(n))
  if(!sd){return(avpc$mm)}else{return(avpc$sd)}
}

avg_pat_diff_clust <- function(data){
  # Data table of metadata
  data %>%dplyr::count(cluster, multi_q, cond_cln) %>% 
    dplyr::group_by(cluster, cond_cln) %>% 
    dplyr::summarise(mp = mean(n), medp = median(n)) %>% 
    dplyr::group_by(cluster) %>%
    dplyr::mutate(avg_diff = diff(mp), med_diff = diff(medp))
}

pull_logfc <- function(perm, feature, col="obs_log2FD"){
  dplyr::filter(
    perm@results$permutation,
    stringr::str_detect(clusters, feature)
  ) %>% dplyr::pull(!!sym(col))
}
