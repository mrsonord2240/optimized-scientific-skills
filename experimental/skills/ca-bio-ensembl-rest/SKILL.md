---
name: ca-bio-ensembl-rest
description: Use when retrieving Ensembl gene, sequence, comparative, or variant data through REST.
category: Data Analysis
license: MIT
author: GPTomics
---

# Ensembl REST

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Gene, sequence, overlap or cross-references | `routes/lookup.md` |
| Variant effects or LD | `routes/variants.md` |
| Orthologs, paralogs or gene trees | `routes/compara.md` |

Record species, assembly and release; live and archive availability may change. Honor service throttling and Retry-After. GPTomics bioSkills, MIT notice in `LICENSE`, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked requests 2.34.2 and Ensembl release 116; current endpoint contracts were not checked here.
