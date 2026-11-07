# Preprocessing problems

| Symptom | Cause | Fix |
|---------|-------|-----|
| HVGs look random, clustering is mush | `seurat_v3` fed log-normalized data, or `seurat` fed raw | Feed each flavor its input; `layer='counts'` for `seurat_v3` |
| A healthy cell type disappeared | Flat mito cutoff removed high-mito parenchyma | MAD with a tissue-aware cap; inspect what was removed |
| Values shrink after rerunning; Scanpy warns X looks log-transformed | Normalized already-normalized data | Normalize raw once; restore from `layers['counts']` |
| `ModuleNotFoundError: skmisc` | `seurat_v3` needs scikit-misc | Install `scikit-misc` |
| QC metrics missing from `.obs` | `calculate_qc_metrics` `inplace` defaults to False | Pass `inplace=True` (the script does) |
| Proliferation or activation signal vanished | Regressed out `total_counts` or cell cycle | Do not reflexively regress |
| New "stressed" cluster | Warm-dissociation IEG and HSP artifact | Score the module; exclude those genes from HVGs |
| Hb or Ig off-target everywhere | Ambient RNA | `routes/ambient.md`, on the raw matrix, before QC |
| Dominant spike deflated other genes | Compositional see-saw | `exclude_highly_expressed=True` or scran |
| Almost everything filtered | MAD is 0 (more than half share a value), so every other cell is an outlier | Stop MAD filtering; use fixed cutoffs when any log QC MAD is 0 or under 80% survive |
| Mito filter removes nothing, `percent.mt` is 0 | `^MT-` pattern against Ensembl IDs | Use symbols or a mito gene list |
| Shallow batch over-filtered, deep batch under-filtered | Global MAD across samples | Per-sample MAD (`batch_key`) |
| Batch effects baked into QC, soup or doublet calls | Merged before steps up to doublet detection | Run them per sample, then merge |

Done when the fix has been applied and the script rerun.
