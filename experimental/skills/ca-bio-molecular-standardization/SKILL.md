---
name: ca-bio-molecular-standardization
description: Use when standardizing molecular structures for registration, deduplication, dataset joining, or QSAR.
category: Data Analysis
license: MIT
author: GPTomics
---

# Molecular Standardization

Choose a route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Need | Read |
|---|---|
| ChEMBL-style representation and parent selection | `routes/chembl.md` |
| Explicit RDKit, tautomer or isotope policy | `routes/rdkit.md` |
| Deduplicated QSAR data and replicate review | `routes/qsar.md` |
| Related analysis | `routes/related.md` |

Choose and record fragment, charge, isotope, stereo and tautomer policies. Retain flagged or excluded input status; ChEMBL parent selection alone does not canonicalize tautomers.
Source checked: RDKit 2024.09+, chembl_structure_pipeline 1.2+.
