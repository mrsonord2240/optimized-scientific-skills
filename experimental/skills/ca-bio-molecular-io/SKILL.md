---
name: ca-bio-molecular-io
description: Use when reading, writing, or converting molecular structure files.
category: Data Analysis
license: MIT
author: GPTomics
---

# Molecular I/O

Choose a route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Need | Read |
|---|---|
| Parse or write SMILES and SDF | `routes/read-write.md` |
| Convert MOL2 or PDBQT | `routes/convert.md` |
| InChI identity and tautomer distinctions | `routes/identity.md` |
| Stereo, aromaticity, bond-order or parse errors | `routes/diagnose.md` |
| Related analysis | `routes/related.md` |

Reject blank and zero-atom molecules, retain failed records, and verify representation and stereochemistry through conversion.
Source checked: RDKit 2024.09+, Open Babel 3.1.1+.
