# Seurat v5 bridge integration

```bash
Rscript scripts/seurat_bridge_integration.R rna.rds multi.rds atac.rds bridge_query.rds
```

Optional trailing arguments: `ndims` (20), `first_lsi_dim` (2), `SCT|LogNormalize` (LogNormalize). If Seurat and Signac live in a private library, set `BIO_SKILLS_R_LIB=/path/to/R-lib` before `Rscript` (PowerShell: `$env:BIO_SKILLS_R_LIB = 'C:\path\to\R-lib'`).

Inputs, all `saveRDS` Seurat objects: `rna.rds`, the labelled scRNA reference (`meta.data$celltype`; NormalizeData, ScaleData, RunPCA, and `RunUMAP(return.model = TRUE)`); `multi.rds`, the paired multiome bridge with `RNA` (normalized) and `ATAC` (RunTFIDF, RunSVD giving `lsi`); `atac.rds`, the unpaired scATAC query on the same peak set as the bridge ATAC assay, RunTFIDF only. Output: the query with `predicted.celltype`, `predicted.celltype.score` and `ref.umap`.

- `first_lsi_dim=2` drops LSI_1: keep it only if DepthCor confirms LSI_1 tracks depth, else pass `1`.
- Pass `SCT` if the reference and bridge RNA were SCTransformed.
- `MapQuery` takes the object `PrepareBridgeReference` returned as `reference`, not the original scRNA object (the script does).
- Check `predicted.celltype.score` before trusting the labels: on synthetic data 81% of query cells recovered their type with LogNormalize and 75% with SCT (chance 33%). Rare query-only populations get mislabelled.

Done when `bridge_query.rds` exists and the answer reports the score distribution.
