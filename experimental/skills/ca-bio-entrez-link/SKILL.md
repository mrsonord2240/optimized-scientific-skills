---
name: ca-bio-entrez-link
description: Use when finding cross-database links between known NCBI Entrez records.
category: Data Analysis
license: MIT
author: GPTomics
---

# Entrez links

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Discover or retrieve links for a few source UIDs | `routes/link.md` |
| Link a large UID list through history | `routes/history.md` |
| Chain links or obtain related-record scores | `routes/chain.md` |

Confirm source UID and database together, choose the linkname deliberately, and preserve link direction. ELink relationships are not guaranteed symmetric. GPTomics bioSkills, MIT, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked Biopython 1.83+ and Entrez Direct 21.0+; live link tables may change.
