---
name: ca-bio-atac-seq-footprinting
description: Use when detecting or comparing transcription-factor footprints in ATAC-seq data.
category: Data Analysis
license: MIT
author: GPTomics
---

# ATAC-seq TF footprinting

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Two-condition motif footprints and differential activity | `routes/tobias.md` |
| Bulk or per-cluster multiscale footprints | `routes/scprinter.md` |
| Independent footprint calls and concordance | `routes/concordance.md` |

Use a matching reference, consensus peaks and blacklist; assess depth and Tn5 bias. A footprint is evidence of protection, and its absence does not prove absence of binding. GPTomics bioSkills, MIT notice in `LICENSE`, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked TOBIAS 0.17.5, HINT-ATAC 1.0.2, pyDNase 0.3.0 and scPrinter 1.2.0.
