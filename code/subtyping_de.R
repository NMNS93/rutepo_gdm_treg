library(Seurat)
library(metap)

subtyping_de <- function(x, label=""){
  # Perform SCT cluster identification and return DE table 
  ident_set <- as.numeric(as.character(unique(x$seurat_clusters)))
  markers <- lapply(
    ident_set,
    function(ident){
      ms <-  FindConservedMarkers(x, assay = "SCT", ident.1=ident, grouping.var = "control",verbose = FALSE)
      ms$ident = ident
      ms$label = label
      ms$gene = rownames(ms)
      return(ms)
    }
  )
  res <- Reduce(rbind, markers)
  return(res)
}


add_cluster_labels <- function(x, cluster_csv){
  labels = readr::read_csv(cluster_csv)
  rownames(labels) = labels$seurat_clusters
  labels$seurat_clusters = factor(labels$seurat_clusters)
  x_md = data.frame(seurat_clusters=x@meta.data$seurat_clusters, cc = Cells(x))
  meta = left_join(x_md, labels)
  rownames(meta) = meta$cc
  x = AddMetaData(x, meta)
  return(x)
}

get_cluster_average <- function(v){
  Idents(v) = "cluster"
  levels(v) = unique(v$cluster)
  avpv = AverageExpression(v, assays="SCT", group.by = "cluster", return.seurat = T)
  avpv$cluster = Cells(avpv)
  Idents(avpv) = "cluster"
  return(avpv)
}

add_label_data <- function(x, label_csv){
  d = readr::read_csv(label_csv)
  rownames(d) = as.character(d$cluster)
  x = AddMetaData(x, d)
  return(x)
}

de_top_n_genes <- function(de, n){
  dplyr::group_by(de, ident) %>% dplyr::top_n(n, abs(CONTROL_avg_log2FC)) %>% pull(gene)
}

filter_for_top_markers <- function(x){
  x %>%
    dplyr::filter(!(grepl(gene, pattern = "RP[SL]"))) %>%
    dplyr::filter(!(grepl(gene, pattern = "MT-"))) %>%
    group_by(label, ident) %>%
    dplyr::top_n(10, -log10(CONTROL_p_val_adj) * -log10(CASE_p_val_adj)) %>%
    dplyr::filter(abs(CONTROL_avg_log2FC) > 0.25 & abs(CASE_avg_log2FC)>0.25)
}

reorder_by_expr <- function(genelist, data, subtype_order, NN=10){
  data = data[genelist,]
  sets = lapply(
    subtype_order,
    function(i)  order(-data[,i])[1:NN]
  )
  gsel = unique(unlist(sets))
  gtop = rownames(data)[gsel]
  gend = rownames(data)[ which(!rownames(data) %in% gtop)]
  genelist2 = c(gtop, gend)
  return(genelist2)
}