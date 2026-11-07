# Seurat marker ranking
~~~sh
Rscript scripts/seurat_markers.R clustered.rds markers.csv
~~~
The source uses FindAllMarkers with positive markers, Wilcoxon, logfc.threshold=0.25 and min.pct=0.1. Review avg_log2FC and pct.1/pct.2; the script writes all returned markers and does not infer labels. Set identities and normalized expression appropriately before running. Check installed Seurat API.
**Done when:** the table and expression specificity are inspected for each cluster.
