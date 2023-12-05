# Use scRepertoire to preprocess VDJ results
library(Seurat)
library(scRepertoire)
library(dplyr)
library(ggplot2)
library(data.table)
library(parallel)
source("code/preprocessing_qc.R")
# 
# # Variables
# outdir = fs::dir_create(here::here("data/seurat/220802_Shangaris"))
# figdir = fs::dir_create(here::here("figures/220802_seurat"))
# meta_csv = here::here("data/metadata/220721_metadata_cln.csv")
# vdj_locs = here::here("data/cellranger/220907_Shangaris_CellRangerDirs_vdj.csv")
# 
# # Data
# tenx = qs::qread(fs::path(outdir, "seurat_tenx_clustered.qs"))
# DefaultAssay(tenx) = "RNA"
# tenx = AddPanicosMeta(tenx, meta_csv, cols_regex="developed|condition|control")
# 
# # There are duplicate cell barcodes present across sequencing libraries. This may
# # be because the same library kit was used or mixed for tagging cells.
# # To overcome this, we process the single cell data one sequencing run at a time.
# contig_by_patient <- function(sc, vdj_csv, run_id){
#   stopifnot(
#     "orig.ident" %in% names(sc@meta.data) &
#       "HTO_maxID" %in% names(sc@meta.data)
#   )
#   # Load VDJ data
#   contig = read.csv(vdj_csv)
#   # Filter input seurat object to cells of interest based on run_id (must be orig.ident field)
#   run_cells = subset(sc, cells=WhichCells(sc, expression = orig.ident == run_id))
#   # Rename cells of interest if any barcode suffixes are present
#   new.cells = stringr::str_remove(Cells(run_cells), "_[0-9]+$")
#   upcells = RenameCells(run_cells, new.names=new.cells)
#   # Create contig list
#   htl = createHTOContigList(contig, upcells, group.by="HTO_maxID")
#   # Optional: Strip prefixes from barcodes with stripPrefix (not required)
#   return(htl)
# }
# 
# # Load VDJ contigs
# vdj_recs = readr::read_csv(vdj_locs, show_col_types = F) %>% split(., .$sample_id)
# contig_list = mclapply(
#   vdj_recs,
#   function(x){
#     contig_by_patient(tenx, x$vdj_csv, x$sample_id)
#   },
#   mc.cores=10
# )
# 
# get_contig_condition <- function(contig_list, tenx){
#   annos = tenx@meta.data %>% dplyr::select(HTO_maxID, condition) %>% unique()
#   df = data.frame(HTO_maxID=factor(names(contig_list), levels=names(contig_list)))
#   df = left_join(df, annos)
#   stopifnot(identical(df$HTO_maxID, names(contig_list)))
#   return(df$condition)
# }
# 
# # Format and add condition annotation
# contig_list = Reduce(c, contig_list)
# contig_condition = get_contig_condition(contig_list, tenx)
# 
# # Combine TCR records for all cells and samples
# cl = combineTCR(contig_list, samples=names(contig_list), ID=contig_condition, cells="T-AB")
# qs::qsave(cl, fs::path(outdir, "screp_combined_tcr.qs"))
# 
# # TODO: Generate plots
# # From here on, addVariable and subsetContig allow us to play with the TCR VDJ data
# quantContig(cl, cloneCall="strict", scale=F) + coord_flip() + theme(legend.position = "none")
# quantContig(cl, cloneCall="strict", scale=T) + coord_flip() + theme(legend.position = "none")
# abundanceContig(cl)
# #abundanceContig(cl) + facet_wrap(stringr::str_split(v$data$values, "_", simplify=T)[,2])
# lengthContig(cl, cloneCall="aa", chain = "both")   # large contigs meaningful?
# compareClonotypes(cl, 
#                   numbers = 2, 
#                   cloneCall="aa", 
#                   graph = "alluvial")
# 
# clonalHomeostasis(combined, cloneCall = "gene", 
#                   cloneTypes = c(Rare = 1e-04, 
#                                  Small = 0.001, 
#                                  Medium = 0.01, 
#                                  Large = 0.1, 
#                                  Hyperexpanded = 1))
