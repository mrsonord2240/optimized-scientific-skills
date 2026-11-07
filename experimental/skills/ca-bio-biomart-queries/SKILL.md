---
name: ca-bio-biomart-queries
description: Use when querying BioMart for gene identifiers, coordinates, annotations, or orthologs in bulk.
category: Data Analysis
license: MIT
author: GPTomics
---

# BioMart queries

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Map a list of identifiers or retrieve GO terms | `routes/ids.md` |
| Coordinates or a version-pinned query | `routes/coordinates.md` |
| Ortholog table | `routes/orthologs.md` |

Confirm dataset, species, build, release and available attributes; BioMart can return multiple rows per input. GPTomics bioSkills, MIT, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked pybiomart 0.2.0, biomaRt 2.62.1 and Ensembl release 116; live schemas may change.
