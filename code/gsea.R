library(fgsea)
library(dplyr)
library(rbioapi)
library(data.table)

gsea_runner_v3 <- function(tbl, gmt_files){

  # Split by clusters
  cls = as.numeric(as.character(tbl$seurat_clusters))
  pc <- split(tbl, cls)
  
  tbl_gsea <- mclapply(
    pc,
    function(e){
      # Log
      i = e$seurat_clusters[[1]]
      log_info("GSEA running for seurat_cluster[{i}] on {length(gmts)} sets.")
      
      # Rank genes for GSEA
      vrank= -log10(e$p_val) * sign(e$avg_log2FC)
      names(vrank) = e$gene
      vrank = vrank[order(vrank)]
      
      # Run GSEA
      res <- gsea_helper_v3(vrank, gmt_files)
      res$seurat_clusters = e$seurat_clusters[[1]]
      res$cluster = e$cluster[[1]]
      res$label = e$label[[1]]
      res$condition = e$condition[[1]]
      return(res)
    },
    mc.cores=3
  ) %>% Reduce(f=bind_rows)
  
  
  return(tbl_gsea)
}

gsea_single_v3 <- function(ranks, pathways){
  
  # Run fGSEA. N.B. Recommended that all expressed genes are considered:
  #   https://github.com/ctlab/fgsea/issues/26#issuecomment-1012061210
  fgseaRes <- tryCatch(
    fgsea(pathways, ranks, maxSize=500),
    error=function(e) data.table()
    
  ) # Min genes in set
  
  return(fgseaRes)
}

gsea_helper_v3 <- function(vrank, gmt_files){
  # Run GSEA over mutitple gmt files
  res <- mclapply(
    gmt_files,
    function(o){
      
      # Load pathway set
      pathways = gmtPathways(o)
      
      # Setup ranks. Previously calculated at -log10(padj) * sign(log2fc)
      ranks = vrank
      
      fgseaRes = gsea_single_v3(ranks, pathways)
      
      if(nrow(fgseaRes) > 0){
        fgseaRes$gmt = fs::path_file(o)
        fgseaRes = dplyr::arrange(fgseaRes, padj)
      }
      
      },
    mc.cores=6
  )
  
  res = data.table::rbindlist(res)
  
  return(res)

}

gsea_runner <- function(tbl, gmt_files, field="avg_log2FC"){
  
  # Split by clusters
  pc <- split(tbl, tbl$cluster)
  
  tbl_gsea <- mclapply(
    names(pc),
    function(e){
      log_info("GSEA running for {e} on {length(gmts)} sets.")
      pc_tbl = pc[[e]] %>% dplyr::filter(condition=="CASE")
      res <- gsea_helper(pc_tbl, gmt_files, field)
      res$dataset = pc_tbl$dataset[[1]]
      res$seurat_clusters = pc_tbl$seurat_clusters[[1]]
      res$cluster = pc_tbl$cluster[[1]]
      return(res)
    },
    mc.cores=3
  ) %>% Reduce(f=bind_rows)
  

   return(tbl_gsea)
}

gsea_helper <- function(tbl, gmt_files, field){
  # Wrapper to run GSEA for multiple gmt files
  res <- mclapply(
    gmt_files,
    function(o) gsea_single(tbl, o, field),
    mc.cores=6
  ) %>% Reduce(f = bind_rows)
  return(res)
}

gsea_single <- function(tbl, gmt, field="avg_log2FC"){
  # Table should have DE results in the typical format. Only requirements are
  #   tbl$gene and tbl[[field]]
  # GMT files contain pathways genesets
  
  pathways = gmtPathways(gmt)
  
  # Setup ranks for avgLog2FC
  ranks = tbl[[field]]
  names(ranks) = tbl$gene
  if(field == "p_val_adj"){
    ranks = -log10(ranks)
  }
  ranks = ranks[order(ranks)]

  # It is  recommended that all expressed genes are considered for fgsea:
  #   https://github.com/ctlab/fgsea/issues/26#issuecomment-1012061210
  fgseaRes <- tryCatch(
    fgsea(pathways, ranks, minSize=3, maxSize=500),
    error=function(e) data.table()
  
    ) # Min genes in set
  
  if(nrow(fgseaRes) > 0){
    fgseaRes$gmt = fs::path_file(gmt)
    fgseaRes$field = field
    fgseaRes = dplyr::arrange(fgseaRes, padj)
  }
  
  return(fgseaRes)
}

run_ora <- function(genes){
  go_ <- rba_enrichr(gene_list = genes, gene_set_library = "GO_.*2021")
  kegg_ <- rba_enrichr(gene_list = genes, gene_set_library = "KEGG_2021_Human")
  wiki_ <- rba_enrichr(gene_list = genes, gene_set_library = "WikiPathway_2021_Human")
  combined<- c(go_, list("KEGG_2021_Human" = kegg_, "WikiPathway_2021_Human" = wiki_))
  result <- lapply(names(combined), FUN=function(x){
   gset = combined[[x]]
   if(nrow(gset)==0){return(gset)}
   gset$dataset = x
   gset = dplyr::filter(gset, Adjusted.P.value <= 0.01)
   return(gset)
  }) %>% rbindlist(., fill=T)#Reduce(f=dplyr::bind_rows)
  if(nrow(result)==0){return(result)}
  # Append gene ratio
  grmatrix = stringr::str_split(result$Overlap, "/", simplify=T)
  result$de_genes =  as.numeric(grmatrix[,1])
  result$gene_set_size =  as.numeric(grmatrix[,2])
  result$gene_ratio = result$de_genes / result$gene_set_size
  return(result)
}

plot_ora <-function(x){
  ggplot(x, aes(y=Term, x=-log10(Adjusted.P.value), fill=dataset)) + geom_col() + theme_bw()
}

plot_ora_v2 <- function(x, high.col){
  # x is output from run_ora
  stopifnot(
    all(
      c("mltp", "Term", "gene_ratio", "gene_set_size") %in% names(x)
    )
  )
  x = dplyr::mutate(x,
                    Term = forcats::fct_reorder(Term, mltp))
  ggplot(x,
         aes(x=mltp, y=Term)
  ) +
    geom_point(
      aes(fill=gene_ratio, size=gene_set_size),
      shape=21
    ) +
    scale_x_continuous(
      limits=c(0, NA)
    ) +
    labs(x="-log10(adj_p_value)") + 
    geom_vline(xintercept = -log10(0.05), linetype="dashed", color="grey") +
    theme_bw() +
    scale_fill_gradient(
      low="lightgray", high=high.col
    )
  
}

plot_ora_comb <- function(x, fscale){
  ggplot(x,
         aes(x=mltp, y=forcats::fct_reorder(Term, mltp))
  ) +
    geom_point(
      aes(fill=dataset, size=gene_ratio),
      shape=21
    ) +
    scale_x_continuous(
      limits=c(0, NA)
    ) +
    labs(x="-log10(adj_p_value)") + 
    geom_vline(xintercept = -log10(0.05), linetype="dashed", color="grey") +
    theme_bw() +
    facet_wrap(~label, scales="free_x") +
    scale_fill_manual(values=fscale)
}

unravel_genes <- function(v){
  # Convert a vector of genes concatenated with ';' to a flat string
  unique(stringr::str_split(paste0(v, collapse=";"), ";")[[1]])
}

get_rec_hall <- function(x) {
  t3 = x %>% dplyr::group_by(seurat_clusters) %>% dplyr::filter(stringr::str_detect(pathway, "HALLMARK")) %>%
    top_n(3, -padj) %>% ungroup() %>% group_by(pathway) 
  t3s = t3 %>%
    dplyr::summarise(count=n(), clusters=paste0(cluster, collapse=","), cfam=paste0(seurat_clusters, collapse=","), NESM=mean(NES)) %>%
    top_n(5, count) %>% arrange(count)
  return(list("top3"=t3, "summary"=t3s))
}


detect_many <- function(string, patterns){
  res <- sapply(
    string,
    FUN = function(string) {
      pat = paste0(patterns, collapse="|")
      any(
        stringr::str_detect(string, pat)
      )
    }
  )
  
  return(res)
}

filter_sig_pathway_genes <- function(data, gml){
  res = data %>% group_by(cluster, gmt) %>% top_n(1, -log10(padj)) %>% top_n(1, abs(NES)) %>%
    dplyr::mutate(gmt_short = stringr::str_split(gmt, "_", simplify=T)[,2]) %>%
    dplyr::mutate(setname = paste0(gmt_short, "_", pathway)) %>%
    select(setname, pathway) %>%
    dplyr::mutate(geneset = lapply(pathway, function(i) gml[[i]]))
  
  htc_set = res$geneset
  names(htc_set) = res$pathway
  htc_set = htc_set[unique(res$pathway)]
  
  return(list(res, htc_set))
}

add_escape_gsea_meta <- function(es, x){
  # es is escape ssgea result
  # x is Treg.CD4 object
  stopifnot(all(rownames(es) == rownames(x@meta.data)))
  esd = data.table(es)
  esd$cells = rownames(es)
  esd = cbind(
    esd,
    dplyr::select(x@meta.data, orig.ident, multi_q, seurat_clusters, label, cluster, clust_group)
  )
  return(esd)
}


filter_na_cells <- function(x, cols){
  na_cells = dplyr::select(x@meta.data, all_of(cols))
  na_counts = rowSums(is.na(na_cells))
  dropcells = names(which(na_counts>0))
  selectcells = Cells(x)[!(Cells(x) %in% dropcells)]
  subset(x, cells = selectcells)
}



make_escape_ridge_split <- function(x, single_set, colors, prefix, figdir){
  mc = colorRampPalette(as.character(colors))(12)
  setpre = stringr::str_split(single_set, "_", simplify=T)[,1]
  replot = ridgeEnrichment(x@meta.data, gene.set = single_set, group = "cluster", facet = "cond_cln", add.rug = TRUE, colors=mc)
  spplot = splitEnrichment(x@meta.data, x.axis = "cluster", split = "cond_cln", gene.set = single_set, colors = colors)
  ggsave(fs::path(figdir, paste0(prefix, "_escape_ridge_", setpre, ".png")), replot, width=8, height=10)
  ggsave(fs::path(figdir, paste0(prefix, "_escape_split_", setpre, ".png")), spplot, width=15, height=8)
}
