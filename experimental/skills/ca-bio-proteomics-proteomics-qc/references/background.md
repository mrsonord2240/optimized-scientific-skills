# Background

Not needed to run a route. Why the Skill is shaped this way.

## QC is a three-level funnel and the matrix sees the failure last

| Level | Question | Faults localized |
|-------|----------|------------------|
| 1 Instrument / raw signal | Is the LC-MS hardware performing? | Column, spray or emitter, calibration |
| 2 Identification / run | Did the run identify peptides correctly? | Digestion, handling PTM artifacts, contamination |
| 3 Experiment / quantitative | Are the numbers reproducible and comparable? | Loading, batch, outliers, missingness, swaps |

Faults start at the instrument or in identification, but the protein matrix shows only the downstream symptom. The most integrative metric (% MS2 identified, ID count) is the first alarm and the least specific.

## Median normalization hides loading problems

Median or quantile normalization forces a summary statistic equal across samples. A sample that loaded 3x low is a visible shift in a raw boxplot; after normalization it is shifted up by a constant, the boxplots line up, and its noisy low-abundance signal is stretched into mid-range and enters the differential test. The same applies to TMT channel loading.

## Tools by level

| Tool | Level | Role |
|------|-------|------|
| PTXQC (R, CRAN) | 2-3 | Multi-metric report over MaxQuant txt or mzTab |
| RawTools, RawBeans, rawrr, rawDiag, QuaMeter | 1 | Thermo raw-file and vendor-independent metrics |
| Skyline + Panorama AutoQC | 1, longitudinal | Levey-Jennings, CUSUM, moving-range |
| MSstatsTMT | 3, TMT | Summarization and channel plots |
| pandas, limma matrix QC | 3 | The scripts in this Skill |

## Failure notes

- A two-peaked mass-error histogram is almost always monoisotopic mis-assignment or co-isolation, not calibration drift; a generous tolerance still identifies the mis-assigned precursors, so the ID rate looks fine.
- Base CV formula on log2 data is compressed by about ln2 times the mean log2 intensity (15-20x at typical MaxQuant intensities).
- Pearson r on raw intensity is dominated by a few abundant proteins.
- On a fully confounded synthetic set, mean absolute log2FC of 104 truly changed proteins fell from 1.55 to 0.00 after batch removal: confounded batch cannot be corrected.
- Synthetic loading test: a sample at 0.41x total signal showed 0.62x median.
- DIA-NN "High precision" mode silently median-normalizes, halving median CV against "High accuracy".
- Unseeded PCA gave six different PC3 p-values on a 20x600 matrix; PC1 and PC2 were stable.
