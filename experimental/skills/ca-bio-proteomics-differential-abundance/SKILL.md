---
name: ca-bio-proteomics-differential-abundance
category: Data Analysis
description: Use when testing which proteins differ in abundance between conditions from mass-spectrometry quantification.
tool_type: mixed
primary_tool: limma
license: MIT
author: GPTomics
---

# Differential Protein Abundance

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data | Read |
|-----------|------|
| Protein table from a search engine (MaxQuant proteinGroups-style TSV) or a log2 matrix, 3 to 10 samples per group | `routes/protein-table.md` |
| Same, plus PSM or peptide counts per protein | `routes/psm-counts.md` |
| Label-free, many missing values, proteins seen in one group only | `routes/heavy-missingness.md` |
| Peptide or precursor table (evidence.txt, DIA-NN, MSstats input) | `routes/peptide-table.md` |
| A minimum fold change must be part of the tested claim | `routes/min-fold-change.md` |
| Python only, no R, more than 10 samples per group | `routes/python.md` |

Rules on every route:

- Never impute missing values (Perseus downshift, MinProb, kNN, zero). Missing means low. Proteins seen in one group only are reported as "undetected in <group>", never as a fold change.
- Raw intensities are not normalized. The scripts stop on an un-normalized matrix; rerun with the normalization option the message names.
- Batch, donor or pairing goes in the model as a `batch` column of `samples.csv`. Never feed `removeBatchEffect` output to a test.
- Take significance from the BH-adjusted column of the results file the script writes. Look proteins up by name, never by row number.
- Report the rows the script dropped (`*_dropped.csv` or `*_undetected.csv`) separately from the significant ones. They are not "not significant".
- One method on one route is the answer. Do not rerun Welch, OLS or a second normalization to compare call counts.
- Never filter on `abs(logFC) > c` and `adj.P.Val < 0.05` together and call it a controlled false discovery rate.

After the run, only if the request needs it:

| Need | Read |
|------|------|
| Fold changes for GSEA ranking or a figure | `routes/fold-change.md` |
| An error, or a result that looks wrong | `routes/errors.md` |

Tested with limma 3.62.2, DEqMS 1.24.0, proDA 1.20.0, msqrob2 1.14.1, QFeatures 1.16.0, MsCoreUtils 1.18.0, MSstats 4.14.2, ashr 2.2.63 (R 4.4.3, Bioconductor 3.20); pandas 3.0.5, scipy 1.18.1, statsmodels 0.15.0. Install: `BiocManager::install(c("limma","DEqMS","proDA","msqrob2","QFeatures","MSstats","ashr"))`.
