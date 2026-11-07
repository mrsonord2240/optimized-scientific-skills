# Seurat HTODemux
~~~sh
Rscript scripts/htodemux_run.R input.rds HTO demuxed.rds
~~~
Input is a Seurat object with a raw HTO assay and matching cell barcodes. The source uses CLR margin 2 and positive.quantile 0.99; review tag distributions, negative and multiplet calls, staining and pooling balance. The script preserves all cells and assignments. Verify Seurat version and API.
**Done when:** output, assay, thresholds and per-tag assignment quality have been reviewed.
