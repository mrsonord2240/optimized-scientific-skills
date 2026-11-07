---
name: ca-bio-atac-seq-co-accessibility
description: Use when inferring cis-regulatory co-accessibility from single-cell chromatin data.
category: Data Analysis
license: MIT
author: GPTomics
---

# Co-accessibility

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Peak-by-cell matrix, Cicero links | `routes/cicero.md` |
| Existing ArchR project or paired multiome | `routes/multiome.md` |
| Validate or visualize peak-pair links | `routes/validation.md` |

Co-accessibility is a statistical association, not a physical contact or proof of enhancer function. Validate candidates with orthogonal evidence. Preprocessing is covered by `bio-atac-seq-single-cell-atac`; contact-informed enhancer-gene methods by `bio-atac-seq-enhancer-gene-linking`. GPTomics bioSkills, MIT notice in `LICENSE`, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked Cicero 1.3.x monocle3 branch, monocle3 1.3+, ArchR 1.0.3 and Signac 1.17; SCENIC+ was statically reviewed only.
