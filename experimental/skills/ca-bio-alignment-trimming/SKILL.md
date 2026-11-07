---
name: ca-bio-alignment-trimming
description: Use when filtering or cleaning a multiple sequence alignment for downstream analysis.
category: Data Analysis
license: MIT
author: GPTomics
---

# Alignment trimming

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Phylogenetic alignment | `routes/tree.md` |
| HMM profile or column mapping | `routes/profile.md` |
| Selection or codon-aware processing | `routes/selection.md` |
| Entropy, contamination, or ambiguous columns | `routes/specialist.md` |

Compare downstream trees before and after filtering; trim per locus before concatenation and rebuild partition coordinates. Preserve original alignment and column mappings. GPTomics bioSkills, MIT, alignment/alignment-trimming, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked ClipKIT 2.14.0, trimAl 1.4.1, BMGE 1.12 and 2.0, MACSE 2.07 and Divvier 1.01+.
