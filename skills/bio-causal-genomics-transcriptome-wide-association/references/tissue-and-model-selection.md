# Tissue Selection, S-MultiXcan vs UTMOST, and Single-Cell TWAS

### Tissue Selection Protocol

Tissue choice drives TWAS power and false-positive rate; selecting tissues by inspecting TWAS hit count is circular. Run all three of the following on the GWAS sumstats (independent of any TWAS run) and pick the primary TWAS tissue from the intersection:

1. **Stratified LDSC tissue prioritization** (Finucane 2018 Nat Genet 50:621): `ldsc.py --h2-cts <sumstats> --ref-ld-chr-cts <annot> --w-ld-chr <weights>` against 200+ tissue-specific gene expression annotations
2. **CELLEX** (Timshel 2020 eLife 9:e55851): single-cell tissue / cell-type prioritization on the same GWAS
3. **MAGMA gene-property analysis** (de Leeuw 2015 PLoS Comput Biol 11:e1004219): cheaper substitute when LDSC unavailable

**Operational rule:** Primary TWAS tissue = the tissue with FDR-significant enrichment in at least two of the three methods. Run secondary tissues in S-MultiXcan for cross-tissue replication. Bonferroni for tissue selection alone: 0.05 / 200 annotations = 2.5e-4.

### S-MultiXcan vs UTMOST

| Method | Use case | Rationale |
|--------|----------|-----------|
| S-MultiXcan (Barbeira 2019) | Standard GTEx v8 analysis | Pre-computed MASHR weights; lower compute barrier; PCA-regularised inter-tissue regression |
| UTMOST (Hu 2019 Nat Genet 51:568) | Custom eQTL panel with cross-tissue retraining | Higher power at genes with shared cross-tissue eQTL architecture; group-lasso enforces sparsity across tissues |

Benchmarks: Hu 2019, Barbeira 2019. Choose by panel availability first; the methods recover overlapping but non-identical gene sets.

### Single-Cell and Cell-Type-Resolved TWAS

Bulk-tissue TWAS averages over cell composition; sc-eQTL TWAS recovers cell-type-specific regulation but at lower per-cell-type power.

| Panel | Reference | Cells / tissue |
|-------|-----------|----------------|
| OneK1K | Yazar 2022 Science 376:eabf3041 | PBMC, ~982 donors, 14 cell types |
| HipSci iPSC-eQTL | Kilpinen 2017 Nature 546:370 | iPSC, ~317 donors |
| BLUEPRINT | Chen 2016 Cell 167:1398 | Monocytes, neutrophils, T cells |

**Tooling state (2026):** No fully pre-built scPrediXcan equivalent to MASHR; sc-eQTL weights are panel-specific. Train custom PredictDB or use TIGAR-V2's Bayesian DPR on the sc-eQTL matrix.

**Operational rule:** Run standard bulk-tissue TWAS first; run sc-eQTL TWAS in the prioritized cell type as a secondary analysis; require concordance between bulk and sc results before nominating a cell-type-specific gene. Upstream sc preprocessing: cross-reference single-cell/preprocessing.
