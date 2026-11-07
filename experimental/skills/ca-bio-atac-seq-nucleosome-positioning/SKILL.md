---
name: ca-bio-atac-seq-nucleosome-positioning
description: Use when assessing nucleosome positioning or fragment periodicity from ATAC-seq data.
category: Data Analysis
license: MIT
author: GPTomics
---

# Nucleosome positioning

Choose and read a route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Fragment classes, V-plot or approximate repeat length | `routes/diagnostic.md` |
| Occupancy and nucleosome calls | `routes/occupancy.md` |
| Differential position between conditions | `routes/differential.md` |

Use paired-end, deduplicated, MAPQ-filtered, mitochondrial-read-removed BAMs and matched genome-build features. Inspect a V-plot first; flat bands do not support positioning calls. GPTomics bioSkills, MIT notice in `LICENSE`, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked NucleoATAC 0.3.4, DANPOS3 3.2.4, ATACseqQC and pysam 0.24.1.
