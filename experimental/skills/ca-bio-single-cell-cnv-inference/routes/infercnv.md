# inferCNV
~~~sh
Rscript scripts/infercnv_run.R counts.tsv annotations.tsv gene_order.tsv Tcell,Myeloid infercnv_out
~~~
Counts are raw genes by cells; annotations are cell and group, and gene order is gene, chromosome, start, end. Use confident matched normal reference groups, a consistent genome build and a single patient/sample. The source default is cutoff 0.1 for droplet data, denoising, i6 HMM and four threads. Adapt cutoff for full-length assays and verify installed API. No arbitrary median cutoff establishes malignancy.
**Done when:** outputs, reference identity, QC and tentative calls are reviewed with orthogonal evidence.
