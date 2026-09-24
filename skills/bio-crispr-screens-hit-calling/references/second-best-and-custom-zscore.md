# Second-Best sgRNA Rule and Custom z-score Calling

Read when filtering a hit list for single-guide outliers, or when neither MAGeCK nor BAGEL2 fits the design.

## Second-Best sgRNA Conservative Rule

**Goal:** Reduce false positives from single outlier sgRNAs by requiring the second-most-extreme guide per gene to also be a hit.

**Approach:** For each gene, sort sgRNAs by LFC; require the second-best LFC to exceed a threshold. Rejects genes that depend on one extreme guide.

```bash
python scripts/second_best_lfc.py mageck.sgrna_summary.txt --direction neg -o second_best.tsv
```

`second_best_lfc()` in that script is also importable. It returns per-gene `second_best_lfc` and a `single_guide` flag: a gene with one sgRNA has no second guide to check, so it gets NaN and `single_guide=True` instead of silently reading as a pass.

**Rule:** A high-confidence hit has second-best LFC also passing the threshold. A guide-of-one hit has only one extreme guide and should be flagged for orthogonal validation. This rule predates JACKS and is implicit in MAGeCK RRA but explicit elsewhere. Genes with `single_guide=True` (fewer than 2 sgRNAs in the library) have no second guide to check by construction -- always send these to orthogonal validation rather than treating a NaN second-best LFC as a pass.

## Custom z-score Hit Calling (when standard tools don't fit)

**Goal:** Compute gene-level z-scores when neither MAGeCK nor BAGEL2 fits the experimental design.

**Approach:** RPM-normalize, compute per-sgRNA log2 fold-changes, aggregate to gene level, derive z-score from the null distribution of non-targeting controls (cleanest) or all genes (assumes <40% changing), apply BH correction.

```bash
python scripts/custom_zscore_hit_calling.py counts.tsv --ctrl T0_1,T0_2 --treat T18_1,T18_2 --ntc-prefix NonTargeting -o zscore_hits.tsv
```

`custom_zscore_hit_calling()` in that script is also importable. The null comes from non-targeting genes when `--ntc-prefix` is given, else from all genes.

