# Single-sample Seurat
~~~bash
Rscript scripts/seurat_workflow.R filtered_feature_bc_matrix results_seurat
~~~
Done when QC, doublet calls, cluster resolutions and markers are reviewed. Output is a new directory; the saved object has clusters, not validated cell-type labels. This filtered-input example does not remove ambient RNA, demultiplex, integrate samples or perform sample-level DE. Its fixed gene/mitochondrial limits are illustrative and tissue dependent; inspect gene namespace before MT pattern matching. See [guidance](../references/guidance.md).
