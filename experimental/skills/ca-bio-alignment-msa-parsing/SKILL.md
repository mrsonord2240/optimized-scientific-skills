---
name: ca-bio-alignment-msa-parsing
description: Use when parsing, filtering, or analyzing the content of a multiple sequence alignment.
category: Data Analysis
license: MIT
author: GPTomics
---

# MSA parsing

Choose and read one route, then run its command before writing your own code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Sequence IDs, composition, gaps, conserved columns or consensus | `routes/content.md` |
| Filter rows or remove gappy columns | `routes/filter.md` |
| Weights, effective depth, or coevolving columns | `routes/weights.md` |
| Annotation, residue coordinates, confidence, or special inputs | `routes/special.md` |

Column and residue coordinates are zero-based. Normalize `.` gaps and letter case before comparisons; retain case for A2M/A3M insert states. The source loader recognizes FASTA and Stockholm; specify other formats when using Biopython directly. Provider: GPTomics bioSkills, MIT, alignment/msa-parsing, commit d91ed3d563019e649dc854c56ccd62551359488a.

Source checked Biopython 1.88, NumPy 2.0.2 and pyhmmer 0.12.3.
