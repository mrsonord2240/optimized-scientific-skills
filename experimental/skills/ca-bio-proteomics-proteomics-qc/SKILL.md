---
name: ca-bio-proteomics-proteomics-qc
category: Data Analysis
description: Use when judging the quality of bottom-up proteomics runs or a protein intensity matrix and deciding which samples to exclude before differential testing.
tool_type: mixed
primary_tool: pandas
license: MIT
author: GPTomics
---

# Proteomics Quality Control

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data or question | Read |
|-----------------------|------|
| MaxQuant `proteinGroups.txt`, not yet normalized: loading failures, contaminants | `routes/raw-signal.md` |
| Protein matrix: replicate correlation, sample swap, CV, missing values | `routes/matrix-metrics.md` |
| Protein matrix: PCA, batch effect, outlier sample | `routes/pca-batch.md` |
| TMT plexes: channel loading balance | `routes/tmt-balance.md` |
| DIA-NN `report.parquet`: RT fit, peak width per run | `routes/diann-level1.md` |
| MaxQuant `txt/` folder, one multi-metric report | `routes/ptxqc.md` |
| Protein count dropped and the cause is unknown, or runs drift over a queue | `routes/diagnose.md` |
| Writing up the exclusions and the QC call | `routes/writeup.md` |

Rules on every route:

- Inspect raw signal before normalizing. Median normalization erases a low load, so start with `routes/raw-signal.md` unless that has been done.
- Un-normalized columns only for the loading check: MaxQuant `Intensity <sample>` (not `LFQ intensity`), DIA-NN `Precursor.Quantity` (not `Precursor.Normalised`, `PG.MaxLFQ`), raw TMT reporter intensities.
- Zeros in search-engine output are missing values, not measurements. The scripts convert them.
- A check that could not run (one sample per group, a singleton batch) is "not measured", never "clean". Say which checks could not run.
- Thresholds seed a per-instrument control chart; none is a universal pass mark. Exclude a sample, never rescale it, and stop and ask before excluding.
- Seed anything stochastic (PCA). The scripts do.
- The differential test is `bio-proteomics-differential-abundance`; normalization and imputation mechanics are `bio-proteomics-quantification`; DIA q-value internals are `bio-proteomics-dia-analysis`. Do not repeat them here.

Tested with pandas 3.0.6, numpy 2.5.3, scipy 1.18.1, scikit-learn 1.9.1, DIA-NN 2.x reports; PTXQC 1.1.5 (not run).
