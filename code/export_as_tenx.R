# Write function to export as 10x
export_as_tenx <- function(seu, outdir="."){
  
  # Setup data
  mtx = seu@assays$RNA@counts
  features = data.table(ENS_ID="ENS_ID", symbol=rownames(mtx))
  barcodes = colnames(mtx)
  
  # Ensure output directory exists
  fs::dir_create(outdir)
  # Setup file names
  mtx_file = path(outdir, glue("matrix.mtx"))
  feat_file = path(outdir, glue("features.tsv.gz"))
  barc_file = path(outdir, glue("barcodes.tsv.gz"))
  
  # Write files
  Matrix::writeMM(mtx, mtx_file)
  data.table::fwrite(features, feat_file, row.names = F, 
                     col.names = F, sep='\t', quote=F)
  readr::write_lines(barcodes, barc_file)
  
  # To read, simply give Seurat::Read10X(outdir)
  
}

export_meta <- function(seu, outdir="."){
  
  # Ensure output directory exists
  fs::dir_create(outdir)
  meta_file = path(outdir, glue("metadata.tsv"))
  
  # Select fields
  meta = dplyr::select(seu@meta.data,
                barcode=cc, source=orig.ident, id=multi_q,
                condition=cond_cln, cluster_id=seurat_clusters,
                cluster_name=cluster, cluster_group=clust_group,
                dataset=label)
  rownames(meta)=NULL
  
  # Write data
  readr::write_tsv(meta, meta_file)
}

create_geo_directory <- function(){
  sid = rbind(
    cd4@meta.data %>% dplyr::select(multi_q, cond_cln, label) %>% unique(),
    treg@meta.data %>% dplyr::select(multi_q, cond_cln, label) %>% unique()
  )
  sid$title = paste0(
    sid$multi_q, ",", sid$cond_cln, ",", sid$label, ",", "scRNA-seq"
  )
  sid$pp1 = paste0(sid$label, "_barcodes.tsv.gz")
  sid$pp2 = paste0(sid$label, "_features.tsv.gz")
  sid$pp3 = paste0(sid$label, "_matrix.mtx")
  rownames(sid) = NULL
  readr::write_tsv(sid, "data/export/geo_input.tsv")
  
  # Extract the R1 and R2 files for each of these samples
  # Also use this opportunity to create a hard-link for the file
  flocs = readr::read_csv('/genomics/projects/B061-PanicosShangaris/nana/data/cellranger/220907_vdj_file_locs.csv')
  fr1 = purrr::pmap(list(a=flocs$sample_id, b=flocs$tenx_dir), function(a,b) flocs_to_fastq(a,b,R1=T))
  fr2 = purrr::pmap(list(a=flocs$sample_id, b=flocs$tenx_dir), function(a,b) flocs_to_fastq(a,b,R1=F))
  
  
}

flocs_to_fastq <- function(sample_id, tenx_dir, R1=T){
  droot = gsub('processed.*', '', tenx_dir)
  stype = gsub('^.*-', '', sample_id)
  rtype = ifelse(R1, "R1", "R2")
  fs::dir_ls(droot, recurse=T, glob=paste0("*",stype,"*",rtype,"*.fastq.gz"))
}
