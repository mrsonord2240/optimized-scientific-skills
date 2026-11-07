# Background

Not needed to run an analysis. This explains choices the routes and scripts already make.

## Three traps behind most lost-data bugs

Orientation is opposite across ecosystems. AnnData is cells x genes with `obs` on rows; Seurat and SingleCellExperiment are genes x cells. So `adata.X` is the transpose of `LayerData(seu)` and `assay(sce)`. Getting the transpose right but the metadata axis wrong is the commonest silent conversion bug.

Sparse storage compounds it. R `dgCMatrix` is CSC; Scanpy usually stores `X` as CSR. Transposing a CSR matrix gives CSC for free, so a correct AnnData to Seurat hop is a logical transpose plus a CSR-to-CSC change. Dense materialization of 500k cells x 30k float32 genes is about 60 GB.

Conversion is lossy by default: `layers`, `obsp`/`varp`, nested `uns` are dropped and categoricals coerced to character or NA. `adata.raw` has its own `var`, and tools disagree on whether to read `X` or `.raw.X`.

The Cell Ranger filtered matrix is cell-called, not ambient-corrected; the claim that it is "decontaminated" is false.

## Converter evidence from the source

- `zellkonverter::readH5AD(reader='R', raw=TRUE)` returned an SCE with empty `altExpNames` on a real h5ad written by anndata 0.13.x (zellkonverter 1.16.0), although the `.raw` group existed; `schard::h5ad2sce(use.raw=TRUE)` recovered it.
- Seurat to h5ad via `as.SingleCellExperiment()` plus `zellkonverter::writeH5AD()` was verified end to end on a Seurat 5.5.0 PBMC object; the h5ad reloaded correctly in Python, and `scale.data` was absent after the round trip.
- anndataR 1.2.0 was the source's first choice for Python-free h5ad/zarr I/O, with R 4.5 required; sceasy is aging.

## Open points for audit

- The installed R is 4.4.3, so anndataR cannot be installed (its dependency `Rarr` also fails to compile); the converter table keeps it as the preferred route only where R 4.5 exists. The staging notes say zellkonverter's default basilisk environment pins anndata 0.10.9 and cannot read files written by anndata 0.13.
- The source claimed `sc.read_10x_h5` has no `var_names` argument; installed Scanpy 1.12.4 was not rechecked.
- The loader scripts use `gex_only=False`; the source's Seurat example loaded the filtered matrix only.
- `scripts/load_10x.R` adds non-GEX blocks as `CreateAssay5Object` assays; this branch is the script author's, not the source's, and is unexecuted.
