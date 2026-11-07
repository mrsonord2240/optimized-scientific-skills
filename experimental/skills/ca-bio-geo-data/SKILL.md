---
name: ca-bio-geo-data
description: Use when finding or retrieving GEO studies, series matrices, or linked sequencing runs.
category: Data Analysis
license: MIT
author: GPTomics
---

# GEO data

Choose and read a route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Search GEO or links from a publication | `routes/search.md` |
| Inspect a Series and its matrix | `routes/series.md` |
| Supplementary files or SRA runs | `routes/raw.md` |

Check Series relations and platforms before combining samples; submitter matrix processing may be missing or unsuitable for the intended analysis. GPTomics bioSkills, MIT, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked Biopython 1.83+, GEOparse 2.0+ and GEOquery 2.70+; service and schema behavior may change.
