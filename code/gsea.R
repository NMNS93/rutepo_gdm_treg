# Run GSEA analysis

# Library ----
library(fgsea)
library(dplyr)
library(rbioapi)
library(data.table)

# Functions ----

gsea_runner <- function(tbl, gmt_files) {
    # Split by clusters
    cls <- as.numeric(as.character(tbl$seurat_clusters))
    pc <- split(tbl, cls)

    tbl_gsea <- mclapply(
        pc,
        function(e) {
            # Log
            i <- e$seurat_clusters[[1]]
            log_info("GSEA running for seurat_cluster[{i}] on {length(gmts)} sets.")

            # Rank genes for GSEA
            vrank <- -log10(e$p_val) * sign(e$avg_log2FC)
            names(vrank) <- e$gene
            vrank <- vrank[order(vrank)]

            # Run GSEA
            res <- gsea_helper(vrank, gmt_files)
            res$seurat_clusters <- e$seurat_clusters[[1]]
            res$cluster <- e$cluster[[1]]
            res$label <- e$label[[1]]
            res$condition <- e$condition[[1]]
            return(res)
        },
        mc.cores = 3
    ) %>% Reduce(f = bind_rows)


    return(tbl_gsea)
}

gsea_single <- function(ranks, pathways) {
    # Run fGSEA. N.B. Recommended that all expressed genes are considered:
    #   https://github.com/ctlab/fgsea/issues/26#issuecomment-1012061210
    fgseaRes <- tryCatch(
        fgsea(pathways, ranks, maxSize = 500),
        error = function(e) data.table()
    ) # Min genes in set

    return(fgseaRes)
}

gsea_helper <- function(vrank, gmt_files) {
    # Run GSEA over mutitple gmt files
    res <- mclapply(
        gmt_files,
        function(o) {
            # Load pathway set
            pathways <- gmtPathways(o)

            # Setup ranks. Previously calculated at -log10(padj) * sign(log2fc)
            ranks <- vrank

            fgseaRes <- gsea_single(ranks, pathways)

            if (nrow(fgseaRes) > 0) {
                fgseaRes$gmt <- fs::path_file(o)
                fgseaRes <- dplyr::arrange(fgseaRes, padj)
            }
        },
        mc.cores = 6
    )

    res <- data.table::rbindlist(res)

    return(res)
}

