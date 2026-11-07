# Background

Not needed to run an analysis. This explains choices the routes and scripts already make.

## Two hidden assumptions

Normalization assumes near-constant total mRNA per cell. Shifted-log, scran and sctransform all divide by a size factor meant to capture capture efficiency and depth only. That is false for plasma cells, secretory epithelia, hepatocytes, large neurons and S/G2M cells, which carry 2-10x more mRNA; dividing them to a common total deflates every gene that is not one of their few dominant transcripts and partly erases proliferation biology.

QC metrics are biology in disguise. `pct_counts_mt` conflates apoptosis, genuine metabolic activity (cardiomyocytes, hepatocytes and muscle are constitutively 20-40% mito and healthy), dissociation stress and contamination. A flat global cutoff deletes whole healthy populations and the survivors still cluster cleanly, so the loss is invisible. A tidy UMAP proves nothing.

## QC thresholds

| Metric | Value | Note |
|--------|-------|------|
| `min_genes` | 200 | Below is mostly empty droplets or debris |
| log total counts, log genes | 5 MAD | Looser than scater's 3 MAD to avoid cutting biology; log scale because depth is right-skewed |
| top-20-genes fraction | 5 MAD | Flags low-complexity or dying cells; platelets, megakaryocytes and erythrocytes are constitutively high |
| `pct_counts_mt` | 3 MAD plus tissue cap | Nuclei and unknown: MAD only; PBMC 8; cardiac and hepatic about 30; skeletal muscle about 40 |
| `min_cells` (genes) | 3 | |

Fixed cutoffs are a fast first pass that silently delete valid populations; miQC (a mixture of mito fraction and detected genes) helps when that relationship varies across samples.

## Ambient tools compared

SoupX: fast, interpretable rho; fails on homogeneous data; one global soup. CellBender: deep generative, also calls cells, strong on nuclei; over-removes at high FPR, black box, slow. DecontX: no raw matrix needed; relies on cluster purity. Subtracting a shared soup vector can manufacture negative correlations and zero out genes cells genuinely lacked.

## Open points for audit

- Installed Scanpy is 1.12.4 and scikit-misc 0.5.2 (source tested 1.10+). `highly_variable_genes(batch_key=...)` with `seurat_v3` and the `sc.experimental` Pearson-residual HVG location were not rechecked against 1.12.
- The source Scanpy example loaded `var_names='gene_ids'` and matched mito on `var['gene_symbols']`; `.h5` input gives symbols as `var_names` instead. The script handles both.
- `sc.pp.highly_variable_genes(flavor='seurat_v3', batch_key=...)` writes per-batch flags; the script copies only the combined `highly_variable` flag.
