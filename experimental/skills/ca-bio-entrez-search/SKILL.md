---
name: ca-bio-entrez-search
description: Use when searching NCBI databases for identifiers, counts, or searchable fields.
category: Evidence Insight
license: MIT
author: GPTomics
---

# Entrez search

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Search one database or retain a history set | `routes/search.md` |
| Inspect database fields or spelling | `routes/discovery.md` |
| Compare counts across databases | `routes/cross-database.md` |

Set a contact email, inspect the translated query and honor service limits. ESearch returns UIDs rather than records; use `bio-entrez-fetch` to retrieve content. GPTomics bioSkills, MIT, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked Biopython 1.83+ and Entrez Direct 21.0+; service behavior may change.
