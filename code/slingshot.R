library(slingshot)
library(tradeSeq)
library(logger)
library(dplyr)
source("code/plots.R")

run_slingshot <- function(seu, start.clus=1, colors=NA){
  # Run slingshot on TSNE embeddings at cluster start
  emb = Embeddings(seu, reduction="tsne")
  log_info("Calling slingshot")
  ss = slingshot(emb, clusterLabels = seu$seurat_clusters, start.clus=start.clus, reducedDim="tsne")
  
  # Calculate pseudotime
  log_info("Calculating pseudotime")
  pt <- slingPseudotime(ss, na=F)
  
  # Dataset from which we can generate a plot of curves and minimum spanning tree
  log_info("Generating MST plotter")
  ts = data.frame(emb)
  ts$sc = seu$seurat_clusters
  ts$scc = seu$cluster
  
  # Return results
  return(list(ss=ss, pt=pt, ts=ts))
}

annotate_pseudotime <- function(seu, pt){
  log_info("Annotating pseudotime")
  pt = data.frame(pt)
  dim(seu@meta.data)
  dim(pt)
  rownames(pt) = rownames(seu@meta.data)
  seu <- AddMetaData(seu, pt)
  return(seu)
}

run_trade_seq <- function(seu, ss, pt){
  log_info("Running tradeSeq")
  # seu: seurat object; ss: slingshot object; pt:pseudotime object
  counts = seu@assays$SCT@counts
  
  # Get sling curve weights
  cellWeights <- slingCurveWeights(ss)
  
  # Select only SCT variable genes and extract counts
  varfeats = rownames(seu@assays$SCT@scale.data)
  counts = seu@assays$SCT@counts[varfeats,]
  
  # Optional: Add batch to the GAM
  
  # Determine knots to use
  log_info("Determining knots")
  icMat <- evaluateK(counts = counts, sds = ss, k = 3:7, nGenes = 100,
                     verbose = T, plot = FALSE)
  knots_sol = names(which.min(colMeans(icMat)))
  knots_sol = stringr::str_remove(knots_sol, "k: ") %>% as.numeric()
  
  # Setup parallel workers
  BPPARAM <- BiocParallel::bpparam()
  BPPARAM$workers <- 10 
  
  # Call tradeseq GAM model
  log_info("Calling GAM")
  tseq <- fitGAM(counts = counts, pseudotime = pt, cellWeights = cellWeights,
                 conditions = factor(seu@meta.data$cond_cln), # Conditions important to fit to seperate trajectories
                 nknots = knots_sol, verbose = T, parallel=TRUE, BPPARAM = BPPARAM)
  
  # Run association test for dynamic assoc
  log_info("Running association test")
  ATres <- associationTest(tseq, lineages = TRUE, l2fc = log2(2))
  
  return(list(tseq=tseq, atres=ATres, icMat=icMat))
}

# Then run differential progression using KS test. Note this can already be run on the dataset using lineage annots 
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html

run_lin_ks_test <- function(seu){
  lins = purrr::keep(names(seu@meta.data), ~stringr::str_detect(.x, "^Lineage"))
  
  # # And we can ks.test all lineages
  resd = lapply(lins, function(i){
    res = ks.test(
      seu@meta.data[seu@meta.data$cond_cln=="CONTROL", i],
      seu@meta.data[seu@meta.data$cond_cln!="CONTROL", i]
    )
    resdf = data.frame(
      lin = i,
      stat=res$statistic,
      pval=res$p.value,
      method=res$method,
      alt=res$alternative
    )
    return(resdf)
  })
  resf = do.call(rbind, resd)
  resf$pval.adj = p.adjust(resf$pval, method="BH")
  # Error: ks.test is over-bearing
  return(resf)
}


# ATres gives us the genes significantly varied from one lineage to another
# gns = ATres %>% filter(!is.na(pvalue_5) & pvalue_5 <=0.05) %>% rownames
