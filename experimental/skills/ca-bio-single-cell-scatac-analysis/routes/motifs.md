# chromVAR deviations
~~~sh
Rscript scripts/run_chromvar.R peaks.rds results/out cluster1 cluster2 seurat_clusters
~~~
The source helper expects a hg38 peaks assay and matching genome, JASPAR2020 and chromVAR dependencies; output directory must exist. It uses GC/accessibility-matched backgrounds and seed 1; it densifies the count matrix and may require substantial memory. A motif points to a TF family, not direct occupancy or a unique active TF. The differential comparison is descriptive unless the design accounts for independent biological replicates. Verify the process exits successfully and outputs are fresh and complete.
**Done when:** deviations, background, metadata, process status and independent validation are reviewed.
