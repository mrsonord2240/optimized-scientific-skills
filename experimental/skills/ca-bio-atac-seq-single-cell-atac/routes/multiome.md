# Multiome and specialist analyses

```bash
Rscript scripts/multiome_wnn.R paired_multiome.rds integrated.rds
```

Input must be a Seurat object with paired RNA and ATAC assays for the same cells. The command runs the source PCA, LSI and WNN workflow and saves an RDS; inspect modality weights. Cell-cycle and sex QC, PEAKVI mapping and cluster-level chromBPNet require additional inputs and source-specific adaptation described in `references/specialized-topics.md`; they are not performed by this entry point. For linkage use `bio-atac-seq-co-accessibility` or `bio-atac-seq-enhancer-gene-linking` as appropriate.

Done when the requested integration or specialist result is produced with matching features, barcodes, genome and limitations recorded.
