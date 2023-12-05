RunhdWGCNAPipeline <- function(x, xname, cinterest){
  log_info("Setting up")
  stopifnot(cinterest %in% x$cluster)
  features <- rownames(GetAssayData(x, slot='scale.data', assay='SCT'))
  obj <- SetupForWGCNA(x, features=features, wgcna_name=paste0(xname, "_", cinterest))
  
  # We want to construct metacells for each sample (multi_q), and cluster ("cluster") separately.
  # See tutorial on metacells for more details: https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html#construct-metacells
  log_info("Constructing metacells by patient and cluster")
  mobj <- MetacellsByGroups(
    obj, group.by=c("multi_q", "cluster"), reduction="tsne",
    k=25, max_shared=12, min_cells=50, ident.group="cluster", slot="scale.data", assay="SCT")
  
  log_info("Subsetting expression matrix by cluster of interest")
  mobj <- SetDatExpr(mobj, group_name=cinterest, group.by="cluster", assay="SCT", slot="data")
  
  # Test different soft powers:
  log_info("Testing SoftPowers")
  mobj <- TestSoftPowers(
    mobj,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
  )
  # Select power parameter
  power_table <- GetPowerTable(mobj)
  psel = power_table %>% dplyr::filter(SFT.R.sq >= 0.8) %>% dplyr::slice(1) %>% pull(Power)
  qs::qsave(mobj, "test_mobj.qs")
  
  # Construct coexpr network with softpower result
  log_info("Constructing Network")
  coexpr <- ConstructNetwork(
    mobj, soft_power=psel,
    setDatExpr=FALSE,
    overwrite_tom=T,
    tom_name = paste0(xname, "_", cinterest) # name of the topoligical overlap matrix written to disk
  )
  
  # Compute module eigengenes in single cells. These summarise the expression matrix of an entire module.
  log_info("Calculating EigenGenes")
  coexpr <- ModuleEigengenes(coexpr)
  qs::qsave(coexpr, "test_coexpr.qs")
  
  # Find pairwise correlations between genes in single cells and our new eigengenes
  log_info("Getting connectivity")
  conn <- ModuleConnectivity(
    coexpr,
    group.by="cluster", group_name=cinterest
  )
  conn <- ResetModuleNames(conn, cinterest)
  
  # Add hMEs to the object for downstream vis
  log_info("Annotating Module eigengenes")
  hME2 <- GetMEs(conn)
  mods <- colnames(hME2); mods <- mods[mods!="grey"]
  conn@meta.data <- cbind(conn@meta.data, hME2)
  
  # Apply a score using UCell (Alternative summary from using the module eigengenes)
  #conn <- ModuleExprScore(
  #  conn,
  #  n_genes = 5,
  #  method='UCell'
  #)
  
  # Return the completed co-expression matrix object
  return(conn)
}