# Multiome RNA + ATAC, WNN

```bash
Rscript scripts/multiome_wnn.R multiome.rds multi_out lsi_start=2
```

`multiome.rds`: Seurat object with `RNA` counts and an `ATAC` ChromatinAssay on one peak set common to all cells. The script runs RNA PCA, ATAC TF-IDF, `FindTopFeatures(min.cutoff='q0')`, SVD (LSI), writes a `DepthCor` plot, runs `FindMultiModalNeighbors`, UMAP and clustering (`algorithm = 3`, SLM, not Leiden), and saves `multi_out.rds`, `multi_out_depthcor.pdf` and `multi_out_weights.csv`. Options: `rna_dims=30`, `lsi_end=30`.

- `lsi_start=2` drops LSI_1 and is right only if the DepthCor plot shows LSI_1 tracking depth; pass `lsi_start=1` otherwise.
- Merging Multiome datasets needs one peak set: re-quantify all cells against unified peaks, or peak-boundary differences create spurious batch structure.
- The ATAC gene-activity matrix is an approximation of RNA, not measured RNA; do not treat it as the RNA modality.
- Report the weights; a dominant modality is a red flag.

Done when `multi_out.rds` exists and the lsi dimensions used are stated.
