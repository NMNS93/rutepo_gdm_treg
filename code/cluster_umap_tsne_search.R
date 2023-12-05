# Run grid search on clusters
library(Seurat)
library(logger)
library(ggplot2)
library(glue)
library(parallel)
library(png)
library(grid)
library(gridExtra)

search_umap <- function(obj, prefix, reclust_dir, cores=8){
  # Warning: Minimu neighbors and components is 2
  neighbors = c(5, 20) # Neighbors must be >2
  components = c(seq(2,30, by=1)) # Components must be >1
  g = expand.grid(neighbors, components)
  recs = seq(nrow(g))
  log_info("Running UMAP test")
  
  # Special fig dir
  umap_fig_pre = fs::path(reclust_dir, prefix)
  
  mclapply(
    recs,
    function(x){
      neighbors = g[x,1]
      components = g[x,2]
      obj <- Seurat::RunUMAP(obj, n.components=components, n.neighbors=neighbors, dims=1:50)
      plt_title = glue("umap_neigh{neighbors}_comp{components}")
      uplt <- DimPlot(obj, reduction="umap") + labs(title=plt_title)
      ggplot2::ggsave(glue("{umap_fig_pre}_{plt_title}.png"), uplt, width=5, height=5)
    },
    mc.cores=cores
  )
  
  log_info("UMAP test complete")
}

search_tsne <- function(obj, prefix, reclust_dir, cores=8){
  
  # Special fig dir
  tsne_fig_pre = fs::path(reclust_dir, prefix)
  log_info("Running TSNE perplexity")
  
  perplexity =seq(1,50,by=5)
  mclapply(perplexity,
           function(x){
             obj <- Seurat::RunTSNE(obj, perplexity=x, dims=1:50)
             plt_title = glue("tsne_perplexity{x}")
             uplt <- DimPlot(obj, reduction="tsne") + labs(title=plt_title)
             ggsave(glue("{tsne_fig_pre}_tsne_perplexity{x}.png"), uplt, width=5, height=4)
           },
           mc.cores=cores
  )
  log_info("TSNE perplexity complete")
}


bind_figs <- function(reclust_dir, prefix){
  # Funciton to combine figures together
  outfile = fs::path(reclust_dir, glue("{prefix}_figbind.png"))
  pngs = sort(fs::dir_ls(reclust_dir, glob=glue("*{prefix}*")))
  pngs = lapply(pngs, function(x){
    img <- as.raster(readPNG(x))
    rasterGrob(img, interpolate = FALSE)
  })
  nrow = ceiling(length(pngs)/5); ncol= 5
  width = 5*ncol; height=5*nrow
  ggsave(outfile, width=width, height=height,
         marrangeGrob(grobs = pngs, nrow=nrow, ncol=5,top=prefix), limitsize=F)
}

search_clusters <- function(obj, prefix, reclust_dir, cores=8){
  resolution = seq(0.1, 0.8, by=0.05)
  log_info("Cluster test beginning")
  # Special fig dir
  cluster_fig_pre = fs::path(reclust_dir, prefix)
  
  mclapply(
    resolution,
    function(x){
      obj <- FindClusters(obj, resolution=x)
      plt_title <- glue("{prefix}_cluster_resolution{x}")
      uplt <- UMAPPlot(obj) + labs(title=plt_title)
      tplt <- TSNEPlot(obj) + labs(title=plt_title)
      ggsave(glue("{cluster_fig_pre}_umap_{plt_title}.png"), uplt, width=5, height=4)
      ggsave(glue("{cluster_fig_pre}_tsne_{plt_title}.png"), tplt, width=5, height=4)
    },
    mc.cores=8
  )
  
  log_info("Cluster test complete")
}



## Legacy
# Functions to create CLUSTER/UMAP/TSNE plots with multiple parameters
# source(here::here("code/preprocessing_qc.R"))
# library(glue)
# library(ggplot2)
# library(parallel)
# library(logger)
# figout = fs::dir_create(here::here("figures/"))
# 
# # Run PCA and basic clustering on object to prep for tests
# pca_prep <- function(obj, prefix){
#   if(!("pca" %in% Reductions(obj))){
#     obj <- RunPCA(obj, npcs=50, verbose =F)
#   }
#   
#   obj <- RunSeuratClustering(obj)
#   
#   ggsave(
#     fs::path(figout, glue("{prefix}_pca_elbow.png")),
#     ElbowPlot(obj, ndims=50)
#   )
#   
#   return(obj)
# }
# 
# # Rudimentary clustering
# cluster_prep <- function(obj){
#   RunSeuratClustering(obj)
# }
# 
# search_umap <- function(obj, prefix){
#   neighbors = c(5, seq(10,100, by=20), 100)
#   components = c(2, 10, 20, 30, 40)
#   g = expand.grid(neighbors, components)
#   recs = seq(nrow(g))
#   log_info("Running UMAP test")
#   
#   # Special fig dir
#   umap_fig = fs::dir_create(fs::path(figout, "220812_UMAP_test"))
#   
#   mclapply(
#     recs,
#     function(x){
#       neighbors = g[x,1]
#       components = g[x,2]
#       obj <- Seurat::RunUMAP(obj, n.components=components, n.neighbors=neighbors, dims=1:30)
#       uplt <- DimPlot(obj, reduction="umap")
#       ggsave(fs::path(umap_fig, glue("{prefix}_umap_neigh{neighbors}_compon{components}.png")), uplt)
#     },
#     mc.cores=8
#   )
#   
#   log_info("UMAP test complete")
# }
# 
# search_tsne <- function(obj, prefix){
#   
#   # Special fig dir
#   tsne_fig = fs::dir_create(fs::path(figout, "220812_TSNE_test"))
#   log_info("Running TSNE perplexity")
#   
#   perplexity = c(5, 15, 30, 50, 70, 100, 200)
#   mclapply(perplexity,
#    function(x){
#      obj <- Seurat::RunTSNE(obj, perplexity=x, dims=1:30)
#      uplt <- DimPlot(obj, reduction="tsne")
#      ggsave(fs::path(tsne_fig, glue("{prefix}_tsne_perplexity{x}.png")), uplt)
#    },
#    mc.cores=8
#   )
# 
#   log_info("TSNE perplexity complete")
# }
# 
# 
# search_clusters <- function(obj, prefix){
#   resolution = seq(0.1, 1, by=0.1)
#   log_info("Cluster test beginning")
#   # Special fig dir
#   cluster_fig = fs::dir_create(fs::path(figout, "220812_CLUSTER_test"))
#   
#   mclapply(
#     resolution,
#      function(x){
#        obj <- FindClusters(obj, resolution=x)
#        uplt <- UMAPPlot(obj) + NoLegend() | TSNEPlot(obj) 
#        ggsave(fs::path(cluster_fig, glue("{prefix}_cluster_resolution{x}.png")), uplt, width=18)
#      },
#     mc.cores=8
#   )
#   
#   log_info("Cluster test complete")
# }
# 
# run_search_pipeline <- function(obj, prefix){
#   log_info("Running search for {prefix}")
#   obj <- pca_prep(obj, prefix)
#   search_umap(obj, prefix)
#   search_tsne(obj, prefix)
#   
#   obj <- RunDimRed(obj)
#   search_clusters(obj, prefix)
#   log_info("Complete for {prefix}")
# }