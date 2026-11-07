# Preprocess with Scanpy

```bash
python scripts/preprocess_scanpy.py filtered_feature_bc_matrix.h5 preprocessed.h5ad tissue=unknown
```

Input is a Cell Ranger `.h5`, a `filtered_feature_bc_matrix/` directory, or an `.h5ad` with raw counts. The script computes QC metrics (mito, ribo, hemoglobin), flags outliers by MAD, drops them, filters genes, stashes counts in `layers['counts']`, normalizes (`normalize_total` then `log1p`), selects HVGs, runs PCA on the HVGs, and writes `preprocessed.h5ad` plus `<out>.qc.tsv` (MADs, thresholds, survivors).

Options (key=value):

| Option | Default | Use |
|--------|---------|-----|
| `tissue` | `unknown` | `nuclei`, `unknown` = MAD only; `pbmc` caps mito at 8%; `cardiac`, `hepatic` 30%; `skeletal_muscle` 40% |
| `mt_prefix` | `MT-` | `mt-` for mouse |
| `batch_key` | none | An `.obs` column of sample ids: MADs are computed per sample and HVGs per batch |
| `n_top_genes` | 2000 | HVG count |
| `flavor` | `seurat_v3` | `seurat_v3` on raw counts (needs scikit-misc), or `seurat` on log-normalized |
| `min_cells` | 3 | Gene filter |

The script stops, rather than filtering badly, when: no gene matches the mito prefix (Ensembl IDs); a QC MAD is 0 (low-variance, tiny or nuclei data: use documented fixed cutoffs); fewer than 80% of barcodes survive; the input is not integer counts.

- The matrix to load is the filtered one for QC. If ambient correction is wanted, do `routes/ambient.md` first from the raw matrix.
- High mito is ambiguous: apoptosis comes with low gene counts; warm-dissociation stress comes with FOS, JUN, JUNB, EGR1 and HSPA1A/B at normal gene counts and makes a fake "stressed" cluster. Score the dissociation module and exclude those genes from HVGs, or flag them.
- Report removal rates, and what was removed per preliminary cluster or known cell type; high top-20-genes fraction is normal for platelets and erythrocytes.
- Doublet detection runs on raw counts per sample (doublet-detection Skill); use `layers['counts']`.
- Single-nucleus data: standard mito thresholds are meaningless (baseline near zero); rely on counts and genes.

Done when `preprocessed.h5ad` and the QC table exist and you have reported cells in and out.
