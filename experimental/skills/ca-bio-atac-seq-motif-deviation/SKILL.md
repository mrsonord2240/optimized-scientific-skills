---
name: ca-bio-atac-seq-motif-deviation
description: Use when finding which transcription factor motifs vary in accessibility across ATAC-seq samples or single cells with chromVAR.
license: MIT
category: Data Analysis
author: GPTomics
---

# Motif Deviation (chromVAR)

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data and task | Read |
|--------------------|------|
| Bulk peak-by-sample counts, 3 or more samples, motif variability and optional condition contrast | `routes/bulk.md` |
| Seurat/Signac scATAC object, per-cell motif z-scores and cluster markers | `routes/signac.md` |
| ArchR project, per-cell deviations and cluster markers | `routes/archr.md` |
| Bulk time course (5 or more points) | `routes/time-course.md` |
| Reading variability and z-scores, comparing tools or studies | `routes/interpret.md` |
| Footprinting, differential accessibility, co-accessibility, TF-to-target regulons | `routes/other-skills.md` |
| An error, NA z-scores, or a result that looks wrong | `routes/errors.md` |

Rules on every route:

- chromVAR needs a full-scale peakset (5000 or more peaks, typically 50k-200k) and the final peaks. Run motif matching after the peakset is final and rerun if peaks change.
- Background sampling is random: set and report the seed. Keep `niterations=50`, `bs=50`; never go below 30 iterations.
- Bulk: 6 or more samples with real biological variation. With 3-5 samples report variability ranking only; one replicated condition makes chromVAR uninformative (use footprinting or differential accessibility).
- Differential motifs come from limma on z-scores; take significance from `adj.P.Val` (limma has no `FDR` column) and make the first factor level the reference.
- Total library depth per sample is required input, not `colSums(counts)`.
- Report the motif database (JASPAR versus CIS-BP); results differ by database. z-scores are tool-specific: compare across tools or studies only after recomputing on the same peakset and motif set.
- Below 5000 peaks, 1500 total reads per sample (bulk) or 500 cells per cluster (single cell), z-scores are unreliable. All-high variability (above 5) with GC-dominated top motifs signals this: enlarge the peakset or aggregate cells.
- When global accessibility differs more than 5x between cell types, run chromVAR per cluster.

Tested on R 4.4.3 / Bioconductor 3.20: chromVAR 1.30.1, motifmatchr 1.30.0, JASPAR2024 0.99.6, TFBSTools 1.44.0, limma 3.62.1, ArchR 1.0.3, Seurat 5.5.1, Signac 1.17.1 (1.16.0 is the last with `RunChromVAR`).

Derived from GPTomics/bioSkills (Domen Jemec), commit d91ed3d563019e649dc854c56ccd62551359488a (`atac-seq/motif-deviation`); MIT notice in `LICENSE`.
