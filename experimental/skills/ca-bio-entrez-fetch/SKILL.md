---
name: ca-bio-entrez-fetch
description: Use when retrieving full records or summaries for known NCBI database identifiers.
category: Data Analysis
license: MIT
author: GPTomics
---

# Entrez fetch

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Sequence or GenBank records | `routes/sequences.md` |
| Metadata, PubMed, taxonomy, SRA or GEO records | `routes/metadata.md` |
| ClinVar or dbSNP records | `routes/variants.md` |
| Large query result set | `routes/history.md` |

Identify the requester to NCBI, use accession versions where applicable, validate response format and honor service limits. GPTomics bioSkills, MIT, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked Biopython 1.83+ and Entrez Direct 21.0+; live service contracts may change.
