## Multiomic Perturb-seq (RNA + ATAC)

**For chromatin readout:** Use 10X Multiome with CRISPRi/a; sgRNA assignment via the same scATAC-seq library. This section assumes already-quantified RNA (genes x cells) and ATAC (peaks x cells) count matrices for the same cells (shared `obs_names`) -- peak calling from raw fragment files is upstream of this Skill (ArchR or Signac in R, or CellRanger ARC).

```bash
python scripts/multiome_differential.py multiome.h5mu --group "GENE_A KO" --control NTC --out peaks.tsv
```

The script reads `mdata['rna'].obs['mixscape_class']` (per-target call, e.g. 'GENE_A KO', from the Mixscape steps above), propagates it to the
ATAC modality via the shared cell index -- not the pooled `mixscape_class_global`, which merges different target genes' KO cells and dilutes any
perturbation-specific chromatin signal -- and tests KO vs control peaks on `normalize_total` + `log1p` counts with Wilcoxon. Do not use TF-IDF
here: it is for embedding/LSI, and its cell-wise reweighting distorted the Wilcoxon null on a planted-signal synthetic dataset.

```python
# Peak-to-gene linking (which differential peak sits near which differential gene) needs a
# genome annotation file (columns: peak, gene, distance, peak_type):
#   muon.atac.tl.add_peak_annotation(mdata, annotation_file)     # muon 0.1.9: it is in .tl, not .pp
#   muon.atac.tl.rank_peaks_groups(mdata, groupby=..., add_peak_type=True, add_distance=True)
# That puts the nearest-gene / peak-type / distance annotations in
# .uns['rank_genes_groups'] keyed by group -- they are NOT columns of
# sc.get.rank_genes_groups_df. Without an annotation file, report differential genes (PyDESeq2,
# Pertpy Unified Framework section above) and differential peaks (peaks.tsv, above) separately,
# as here.
```
