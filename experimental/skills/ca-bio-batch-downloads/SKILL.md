---
name: ca-bio-batch-downloads
description: Use when downloading batches of records from NCBI with retry or resume requirements.
category: Data Analysis
license: MIT
author: GPTomics
---

# NCBI batch downloads

Choose and read one route, then run its command before writing code. Paths are relative to this skill.

| Need | Route |
|---|---|
| Query results with checkpointed retrieval | `routes/history.md` |
| Known accession or UID list | `routes/ids.md` |
| Bulk genome or gene packages | `routes/datasets.md` |

Identify the requester to NCBI and obey service limits. Compare fetched counts and checksums when available. GPTomics bioSkills, MIT, upstream d91ed3d563019e649dc854c56ccd62551359488a.

Source checked Biopython 1.88, NCBI Datasets CLI 18.37.0 and Entrez Direct 26.0; live limits require current confirmation.
