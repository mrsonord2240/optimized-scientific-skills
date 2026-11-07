---
name: ca-bio-molecular-descriptors
description: Use when calculating molecular fingerprints, physicochemical properties, or shape descriptors.
category: Data Analysis
license: MIT
author: GPTomics
---

# Molecular Descriptors

Choose a route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Need | Read |
|---|---|
| Morgan or MACCS fingerprints | `routes/fingerprints.md` |
| Physicochemical or drug-likeness descriptors | `routes/properties.md` |
| Conformer-dependent 3D descriptors | `routes/shape.md` |
| Partial charges | `routes/charges.md` |
| MAP4 or MHFP6 | `routes/minhash.md` |
| Related analysis | `routes/related.md` |

Standardize inputs and record representation, stereochemistry and parameters. Do not interpret invalid structures as valid all-zero fingerprints.
Source checked: RDKit 2024.09+, numpy 1.26+, pandas 2.2+.
