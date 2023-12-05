# Library to preprocess Shangaris TenX dataset with standard QC parameters

# Dependencies ----
# Filepaths
library(here)
library(fs)
library(logger)
library(qs)

# Data manipulations
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# Single Cell Analysis
library(Seurat)
library(sctransform)
#library(scProportionTest)

# Figs
library(ggplot2)

# Pipelines ----

## Run GDM Preproccing QC Pipeline ----

RunGDMPreproccessingQC <- function(cellranger_csv, hk_csv, meta_csv, cache_dir){

  # Load and integrate data
  log_info("loading TENX objects")
  sobjs <- read_tenx_csv(cellranger_csv)

  log_info("Normalising")
  sobjs <- lapply(sobjs, normalizer)

  log_info("Preprocessing for data integration")
  qs::qsave(sobjs, fs::path(cache_dir, "pre_qc.qs"))
  sobjs <- lapply(sobjs, RunQCPipeline, hk_ref=hk_csv)
  qs::qsave(sobjs, fs::path(cache_dir, "post_qc.qs"))
  
  log_info("Selecting Tregs and CD4s")
  treg <- filter_by_ctype(sobjs, "Treg")
  treg <- lapply(treg, AddPanicosMeta, meta_csv)
  treg <- lapply(treg, applyGDMfilter)
  qs::qsave(treg, fs::path(cache_dir, "treg_preint.qs"))
  
  cd4 <- filter_by_ctype(sobjs, "CD4")
  cd4 <- lapply(cd4, AddPanicosMeta, meta_csv)
  cd4 <- lapply(cd4, applyGDMfilter)
  qs::qsave(cd4, fs::path(cache_dir, "cd4_preint.qs"))

  log_info("Integrating Tregs")
  treg <- RunSCTransformIntegration(treg)
  treg <- QuietTCRGenes(treg)
  DefaultAssay(treg) = "integrated"
  
  log_info("Integrating CD4")
  cd4 <- RunSCTransformIntegration(cd4)
  cd4 <- QuietTCRGenes(cd4)
  DefaultAssay(cd4) = "integrated"
  
  # Cluster and run dimension reduction
  log_info("Running dimensionality reduction")
  treg = RunSeuratClustering(
    RunDimRed(treg, components=2, neighbors=50, perplexity=50),
    resolution=0.35
  )

  cd4 = RunSeuratClustering(
    RunDimRed(cd4, components=2, neighbors=50, perplexity=50),
    resolution=0.2
  )
  
  # Return datasets
  return(list("treg"=treg, "cd4"=cd4))
}

applyGDMfilter <- function(x){
  gdm_control_bool = stringr::str_detect(x$condition, "GDM|CONTROL")
  subcells = Cells(x)[gdm_control_bool]
  x = subset(x, cells=subcells)
  x$cond_cln = ifelse(grepl("GDM", x$condition), "GDM", "CONTROL")
  return(x)
}

# Run TENX QC Pipeline
# Note: QC thresholds decided by exploratory analysis
#' @param obj A Seurat object to filter
#' @param hk_ref A reference file listing housekeeping genes to filter on
RunQCPipeline <- function(obj, hk_ref){

  # Add annotations required for QC filters
  obj = qc_annotate_tenx(obj, hk_ref)
  
  # Remove Doublets and Negative cells using hashtag assignments
  obj = qc_singlets(obj)
  
  # Remove cells 4X away from MAD
  obj = qc_mad(obj)
  
  # Filter cells with <350 expressed genes
  obj = qc_expressed_genes(obj)
  
  # Filter cells with few housekeeping genes expressed
  obj = qc_housekeeping_genes(obj)
  
  # Filter cells with too many mitochondrial genes expressed
  obj = qc_mitochondrial(obj)
  
  # Filter non-T-cell populations
  obj = qc_singleR_tcell(obj)
  
  # Filter runs if required (e.g. duplicates)
  obj = qc_filter_runs(obj)
  
  # Annotate cell cycle scores for downstream regression
  obj = qc_annotate_cell_cycle(obj)
  
  
  return(obj)
}

# Integrate TENX data pipeline using ScTransform v2
## Source: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#tldr
RunSCTransformIntegration <- function(x){
  # Run SCT
  x <- lapply(
    x,
    function(i){
      i = SCTransform(i, vst.flavor = "v2", verbose = FALSE, vars.to.regress=c("S.Score", "G2M.Score")) %>%
        RunPCA(npcs = 30, verbose = FALSE)
    }
  )
  # Prepare
  features <- SelectIntegrationFeatures(object.list = x, nfeatures = 3000)
  x <- PrepSCTIntegration(object.list = x, anchor.features = features)
  # Integrate
  anchors <- FindIntegrationAnchors(object.list = x, normalization.method = "SCT",
                                           anchor.features = features)
  combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  return(combined)
}


# Integrate TENX data pipeline
RunSeuratIntegration <- function(obj_list){
  features <- SelectIntegrationFeatures(object.list = obj_list)
  anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features)
  integrated = IntegrateData(anchorset=anchors)
  DefaultAssay(integrated) = "integrated"
  integrated <- ScaleData(integrated, verbose=F)
  # Optional. split.by argument to scale per dataset as counts were normalised.
  # See Seurat github issues.
  return(integrated)
}

# Functions ----

## Load data -----

#' Read tenx data for a single sample
#' @param proj 10X sequencing/sample run name e.g. L193-B-CD4
#' @param dir 10X project directory containing features, barcodes and matrix files
read_tenx_with_hashtag <- function(proj, dir){
  # Load tenx data
  obj = Read10X(dir)
  sr = CreateSeuratObject(obj$`Gene Expression`, project=proj, assay = "RNA")

  # Create hashtag object
  hto_raw = obj$`Antibody Capture`
  hashtags = stringr::str_extract(rownames(hto_raw), "AHH[0-9]+")
  patient_id = paste0(proj, "-", hashtags)
  rownames(hto_raw) = patient_id
  hto_assay = CreateAssayObject(hto_raw)
  sr[["HTO"]] <- hto_assay
  return(sr)
}

#' Read tenx data into Seurat from project a CSV file
#' @param tenx_locs_csv A file with two header fields: sample_id and tenx_dir,
#'    where sample_id is the project name and tenx_dir is the tenx directory.
read_tenx_csv <- function(tenx_locs_csv, merge=F){
  
  tenx_locs = readr::read_csv(tenx_locs_csv, show_col_types=F)
  proj_c = tenx_locs$sample_id
  dir_c = tenx_locs$tenx_dir
  seurat_obj_list <- purrr::map2(proj_c, dir_c, read_tenx_with_hashtag)
  n_max_obj <- length(seurat_obj_list)
  
  if(merge){
    # Legacy. Should only be used if technical replicates are present with little batch effect,
    # otherwise, use Seurat Data Integration instead.
    # See: https://github.com/satijalab/seurat/issues/1787
    seurat_merged <- merge(
      seurat_obj_list[[1]],
      y=seurat_obj_list[2:n_max_obj],
      add.cell.ids=proj_c,
      project="Shangaris"
    )
    return(seurat_merged)
  }

  return(seurat_obj_list)
}


## QC -------- 

#' Perform pre-integration QC on each seurat object individually
#' @param objects A list of seurat objects
pre_integration_qc <- function(objects, cores=5){
  obs = mclapply(objects, RunQCPipeline, mc.cores=cores)
  return(obs)
}

#' QC annotations
#' 
qc_annotate_cell_cycle <- function(obj){
  CellCycleScoring(obj,
                   s.features = cc.genes$s.genes,
                   g2m.features = cc.genes$g2m.genes)
}

#' Apply various annotations required for QC
qc_annotate_tenx <- function(tenx, hk_ref){
  tenx$ctype <- ifelse(stringr::str_detect(tenx$orig.ident, "CD4"), "CD4", "Treg")
  tenx <- PercentageFeatureSet(tenx, "^MT-", col.name = "percent_mito")
  tenx <- PercentageFeatureSet(tenx, "^RP[SL]", col.name = "percent_ribo")
  tenx <- PercentageFeatureSet(tenx, "^HB[^(P)]", col.name = "percent_hb")
  if(!is.null(hk_ref)){
    parse_hk_genes <- function(hk_csv){
      hk_genes = readr::read_delim(hk_csv, delim=";", show_col_types = F)
      hk_genes = unique(hk_genes$Gene.name)
      return(hk_genes)
    }
    hk_data = parse_hk_genes(hk_ref)
    hk_present <- (intersect(rownames(tenx), hk_data))
    tenx <- PercentageFeatureSet(tenx, features=hk_present, col.name="percent_hk")
  }
  return(tenx)
}

# Filter to singlets based on Hashtag assignment.
qc_singlets <- function(obj){
  obj <- HTODemux(obj)
  obj <- subset(
    obj,
    cells=WhichCells(obj, expression=HTO_classification.global=="Singlet")
  )
  return(obj)
}

# Filter cells with 4X median absolute deviation
qc_mad <- function(obj, med_x=4){
  # Calculate MAD
  med_RNA = median(obj$nCount_RNA)
  mad_RNA_x = stats::mad(obj$nCount_RNA) * med_x
  # Set MAD lower bound
  mad_l = med_RNA - mad_RNA_x
  mad_l = ifelse(mad_l<0, 0, mad_l)
  # Set MAD upper bound
  mad_u = med_RNA + mad_RNA_x

  # Filter cells with outlier counts
  cells_s = WhichCells(obj, expression = nCount_RNA > mad_l & nCount_RNA < mad_u)
  obj = subset(obj, cells=cells_s)
  return(obj)
}

# Filter cells with few expressed genes based on threshold
qc_expressed_genes <- function(x, cutoff=350){
  cells_s = WhichCells(x, expression = nFeature_RNA >= cutoff )
  x = subset(x, cells=cells_s)
  return(x)
}

# Filter cells expressing too few housekeeping genes
qc_housekeeping_genes <- function(x, cutoff=20){
  # Annotate housekeeping genes
  cells_s = WhichCells(x, expression = percent_hk > cutoff)
  x = subset(x, cells=cells_s)
  return(x)
}

# Filter cells with high mitrochondiral genes
qc_mitochondrial <- function(tenx, cutoff=5){
  cells_s = WhichCells(tenx, expression = percent_mito < cutoff)
  tenx = subset(tenx, cells=cells_s)
  return(tenx)
}

# Filter cells that do not have T cell label based on assignments
qc_singleR_tcell <- function(x){
    library(SingleR)
    library(celldex)
    ref <- BlueprintEncodeData()
    scexp <- Seurat::as.SingleCellExperiment(x)
    pred <- SingleR(test=scexp, ref=ref, labels=ref$label.main)
    
    # Assign predictions
    x[["singleR.label"]] = pred$labels

    # Filter any cells that are not identified as T-cells
    x = subset(x, cells=Cells(x)[grepl("T-cells", x$singleR.label)])
    
    return(x)
  }

# Filter runs that are not desired
qc_filter_runs <- function(x, runs=c("L193-B-CD4-AHH05", "L193-B-Treg-AHH05")){
  bool = Reduce(`&`,
                lapply(runs, function(y) x$hash.ID != y)
  )
  subset(x, cells=Cells(x)[bool])
}

## Processing -----

normalizer <- function(x){
    x <- Seurat::NormalizeData(x, assay="RNA")
    # Note: HTO normalisation is included as required downstream for HTODemux
    x <- Seurat::NormalizeData(x, assay="HTO", margin=2)
    x <- Seurat::FindVariableFeatures(x)
    return(x)
}

AddPanicosMeta <- function(tenx, meta_csv, cols_regex = NULL, add_multiq=T){
  if(add_multiq){
  # Setup fields
    tenx$multi_q = stringr::str_replace(tenx$hash.ID, "-CD4|-Treg", "") 
    tenx$run_id = stringr::str_replace(tenx$multi_q, "-AHH.*", "")
  }
  # Load metadata. 
  meta = readr::read_csv(meta_csv, show_col_types = FALSE)
  if(!is.null(cols_regex)){
    meta = dplyr::select(meta, multi_q, matches(cols_regex))
  }
  
  # Add metadata to seurat object
  md = dplyr::left_join(tenx@meta.data, meta, by="multi_q")
  rownames(md) = rownames(tenx@meta.data)  
  tenx <- AddMetaData(tenx, md)
  
  return(tenx)
}

AddSubtypeLabel <- function(obj, sub_name){
  # Label subtypes
  subs = list(
  trsub = tibble(
    seurat_clusters=as.character(seq(0,11)), # TODO: optional input?
    subtype=c("rTreg","pTreg","earlyTCR", "aTreg", "LIMS1_Treg", "Naive Treg", "Th-like", "Unknown", "cytox_Treg", "TRBV7", "rTreg_TRBV30", "Humanin_Treg")
  ),
  cd4sub = tibble(
    seurat_clusters=as.character(seq(0,14)), # TODO: Optional input from file?
    subtype=c("Naive T", "Early Act", "TfH", "Mature T", "Cyto_Interm", "Th17", "T-expansion", "memory Th2", "rT", "Cyto_Late", "Treg", "Naive T2", "Unknown", "Thymocytes_TRBV30", "TRVB7")
  )
  )
  
  add_meta_generic <- function(x, df){
    md = dplyr::left_join(x@meta.data, df)
    rownames(md) = rownames(x@meta.data)  
    x <- AddMetaData(x, md)
    return(x)
  }
  
  obj <- add_meta_generic(obj, subs[[sub_name]])
  return(obj)
}

filter_by_ctype <- function(objs, ctype_filter){
  objs = lapply( objs, function(x){
    tryCatch(
    {
      subset(x, subset = ctype == ctype_filter)
    },
    error=function(cond){return(NA)}
    )}
  )
  objs = objs[!is.na(objs)]
  return(objs)
}

MakeQCFigures <- function(treg, cd4){
  p_featrna = RidgePlot(cd4, "nFeature_RNA", group.by="multi_q") + labs(title="CD4_nFeatureRNA") + NoLegend() |
    RidgePlot(treg, "nFeature_RNA", group.by="multi_q") + NoLegend() + theme(axis.text.y=element_blank()) + labs(y="", title="TREG_nFeatureRNA")
  
  d_ncountrna = rbind(
    dplyr::select(cd4@meta.data, multi_q, nCount_RNA, ctype, control),
    dplyr::select(treg@meta.data, multi_q, nCount_RNA, ctype, control)
  )
  p_d_ncountrna = ggplot(d_ncountrna, aes(x=nCount_RNA)) + geom_histogram() + facet_wrap(control ~ ctype)
  
  no_cells = dplyr::count(d_ncountrna, multi_q, ctype) %>% dplyr::mutate(multi_a=stringr::str_split(multi_q, "-AHH", simplify=T)[,1])
  p_no_cells = ggplot(no_cells, 
         aes(x=n, y=multi_q, fill=multi_a)) + geom_col() + facet_wrap(~ctype) + labs(y="", x="Num_cells", fill="batch")
  
  cond_count = unique(dplyr::select(cd4@meta.data, multi_q, condition)) %>% dplyr::count(condition)
  p_pat_pie = ggplot(cond_count, aes(x="", y=n, fill=condition)) + geom_bar(width=1, stat="identity") + coord_polar("y") +
    scale_y_continuous(n.breaks=11) + theme_minimal() + labs(x="", y="")
  
  qc_plots = list(
    "p_featrna" = p_featrna,
    "p_d_ncountrna" = p_d_ncountrna,
    "p_no_cells" = p_no_cells,
    "p_pat_pie" = p_pat_pie
  )
  
  return(qc_plots)
}

SaveQCFigures <- function(qc_plots, figdir){
  for(i in names(qc_plots)){
    figname = fs::path(figdir, paste0(i, ".png", collapse=""))
    width = 10; height=8;
    #optional: widths and heights conditional on i
    ggsave(figname, qc_plots[[i]], width=10, height=7.5)
  }
}

get_patient_average <- function(v){
  avpv = AverageExpression(v, assays="RNA", group.by = "multi_q", return.seurat = T)
  avpv$multi_q = Cells(avpv)
  avptmeta = unique(v@meta.data[, c("multi_q", "control")])
  rownames(avptmeta) = avptmeta$multi_q
  avpv = AddMetaData(avpv, avptmeta)
  avpv = FindVariableFeatures(ScaleData(NormalizeData(avpv))) # Format
  return(avpv)
}

get_cluster_average <- function(v){
  Idents(v) = "seurat_clusters"
  levels(v) = unique(v$seurat_clusters)
  avpv = AverageExpression(v, assays="RNA", group.by = "seurat_clusters", return.seurat = T)
  avpv$seurat_clusters = Cells(avpv)
  # avptmeta = unique(v@meta.data[, c("seurat_clusters", "cluster_labels")])
  # rownames(avptmeta) = avptmeta$seurat_clusters
  # avpv = AddMetaData(avpv, avptmeta)
  avpv = FindVariableFeatures(ScaleData(NormalizeData(avpv))) # Format
  Idents(avpv) = "seurat_clusters"
  return(avpv)
}


## Clustering -----
# Quiet TCR genes 
# Tim requested this so that clusters aren't driven by TCR status:
# Also previously implemented in another scRepetoire author's package:
QuietTCRGenes <- function(tenx){
  DefaultAssay(tenx) = "integrated"
  tenx <- FindVariableFeatures(tenx)
  vfs = VariableFeatures(tenx)[
    stringr::str_detect(VariableFeatures(tenx), "TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*",
                        negate=T)
  ]
  VariableFeatures(tenx) = vfs
  return(tenx)
}

# Run the standard workflow for Seurat
RunSeuratClustering <- function(obj, resolution=0.3){
  if(!("pca" %in% Reductions(obj))){
    obj <- RunPCA(obj, npcs=50, verbose =F)
  }
  # Updates the seurat_clusters metadata field
  obj <- FindNeighbors(obj, reduction="pca", dims=1:50)
  obj <- FindClusters(obj, resolution=resolution)
  return(obj)
}

# Run dimensionality reduction
RunDimRed <- function(obj, components=20, neighbors=50, perplexity=30){
  if(!("pca" %in% Reductions(obj))){
    obj <- RunPCA(obj, npcs=50, verbose =F)
  }
  obj <- RunUMAP(obj, reduction="pca", dims=1:30,
                 n.components = components, n.neighbors = neighbors)
  obj <- RunTSNE(obj, reduction="pca", dims=1:30, perplexity=perplexity)
  return(obj)
}
## Summaries ----

# 
# n_cells <- function(x) length(Cells(x))
# med_umi_pc <- function(x) median(x@meta.data$nCount_RNA)
# med_genes_pc <- function(x) median(x@meta.data$nFeature_RNA)
# total_genes <- function(x) nrow(x)
# c_sstats <- function(tenx, sstats, label=""){
#   # Calculate summary stats
#   rbind(sstats,
#         data.frame(
#           n_cells=n_cells(tenx),
#           med_umi_pc=med_umi_pc(tenx),
#           med_genes_pc=med_genes_pc(tenx),
#           total_genes=total_genes(tenx),
#           label=label
#         ))
# }

filter_cells_pat_median <- function(x, cutoff){
  cp <- cells_per_pat_clust(x, "x")
  clusts = dplyr::filter(cp, med>=cutoff) %>% pull(seurat_clusters)
  subset(x, seurat_clusters %in% clusts)
}

# Average number of cells per pat per clust.
cells_per_pat_clust <- function(x, label){
 dplyr::count(x@meta.data, multi_q, seurat_clusters) %>% 
  group_by(seurat_clusters) %>% dplyr::summarise(avg=mean(n), sd=sd(n), med=median(n), iqr=stats::IQR(n)) %>% mutate(label=label)
}
