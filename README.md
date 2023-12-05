# RUTEPO scRNA-seq GDM

This repository contains scripts for the analysis of single cell RNA seq data from GDM patients for the Regulatory T Cells in Pregnancy adverse Outcomes (RUTEPO) project.

* `scripts/` - R scripts for analysis numbered in order of processing.
* `code/` - Functions sourced by R scripts to perform analysis.

The pipeline uses Seurat throughout. After quality control, clustering and labelling, we perform proportion tests, differential expression analysis and gene set enrichment. Logsitic regression models are built in the final stage of the pipeline before figures and metrics are calculated. 