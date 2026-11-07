---
name: ca-bio-alignment-msa-statistics
description: Use when calculating identity, conservation, or other statistics from a multiple sequence alignment.
category: Data Analysis
license: MIT
author: GPTomics
---

# MSA statistics

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Pairwise identity matrix | `routes/identity.md` |
| Conservation, entropy, information content, gap statistics | `routes/columns.md` |
| Substitutions, PSSM or alignment quality scores | `routes/scores.md` |
| Corrected distances or alignment uncertainty | `routes/distance.md` |

Read with AlignIO, normalize gap glyphs and letter case before statistics, and inspect alphabet exclusions. Source provider: GPTomics bioSkills, MIT, alignment/msa-statistics, commit d91ed3d563019e649dc854c56ccd62551359488a.

Source checked Biopython 1.88 and NumPy 2.0.2.
