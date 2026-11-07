# Signac LSI from a prepared object
~~~sh
Rscript scripts/lsi_run.R qc_filtered.rds lsi.rds
~~~
Input is a QC-filtered Seurat object with peaks assay and nCount_peaks. Source uses TF-IDF, all features and SVD; the wrapper checks LSI component 1 against depth and excludes it only when absolute correlation exceeds 0.5, then clusters at resolution 0.5. Inspect all component-depth correlations, batch, QC and alternative choices rather than treating this heuristic as universal. Match fragments, peaks and genome build; check package APIs.
**Done when:** output, retained dimensions, depth diagnostics and clusters are reviewed.
