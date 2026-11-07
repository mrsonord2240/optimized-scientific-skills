---
name: ca-bio-atac-seq-differential-accessibility
description: Use when comparing chromatin accessibility between conditions in ATAC-seq data.
category: Data Analysis
license: MIT
author: GPTomics
---

# Differential accessibility

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Replicate BAMs and peak calls, consensus-based comparison | `routes/diffbind.md` |
| Hidden batch or external normalization | `routes/normalization.md` |
| Window, existing count matrix or time-course methods | `routes/alternatives.md` |

Use biological replicates and sample metadata matched to BAM and peak files. Set the case and control explicitly, record normalization and consensus selection, and assess confounding. The preserved script does not invoke a blacklist filter; use prefiltered peaks. GPTomics bioSkills, MIT notice in `LICENSE`, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked DiffBind 3.12+, DESeq2 1.42+, edgeR 4.0+, csaw 1.36+ and related packages.
