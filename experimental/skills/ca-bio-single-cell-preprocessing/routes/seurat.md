# Preprocess with Seurat

```bash
Rscript scripts/preprocess_seurat.R filtered_feature_bc_matrix/ preprocessed.rds tissue=unknown normalization=sct
```

Input is a Cell Ranger directory (`Read10X`) or `.h5` (`Read10X_h5`). The script builds the object (`min.cells=3`, `min.features=200`), computes `percent.mt`, flags MAD outliers (5 MAD on log counts and log genes, 3 MAD on mito, plus the tissue cap), drops them, normalizes, finds variable features and writes `preprocessed.rds`.

| Option | Default | Use |
|--------|---------|-----|
| `tissue` | `unknown` | `nuclei`, `unknown` = MAD only; `pbmc` 8; `cardiac`, `hepatic` 30; `skeletal_muscle` 40 |
| `mt_pattern` | `^MT-` | `^mt-` for mouse |
| `normalization` | `sct` | `sct` runs SCTransform; `lognorm` runs `NormalizeData` (scale factor 10000) then `FindVariableFeatures` |

It stops when `percent.mt` is all zero (symbol pattern against Ensembl IDs), when a MAD is 0, or when fewer than 80% of cells survive.

- Do not regress out `percent.mt`; it is confounded with cell state.
- sctransform v2 is the current form; v1 overfits theta.
- scran deconvolution (R only; needs `quickCluster`; size factors can go negative) suits shallow, high-dropout data; see `routes/normalization.md`.

Done when `preprocessed.rds` exists and you have reported cells in and out.
