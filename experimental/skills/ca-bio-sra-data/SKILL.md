---
name: ca-bio-sra-data
description: Use when resolving public SRA accessions and managing retrieval and validation of archived read files.
category: Data Analysis
license: MIT
author: GPTomics
---
# SRA data access

| Intent | Route |
|---|---|
| Resolve project or study to run IDs | [Metadata](routes/metadata.md) |
| Retrieve public files via ENA | [ENA](routes/ena.md) |
| Retrieve and validate through SRA Toolkit or cloud | [Toolkit](routes/toolkit.md) |

Use only known public data. A run accession does not establish access status; stop on controlled-access uncertainty or authorization failure. See [access and storage notes](references/access.md) and [citations](references/citations.md). This draft was not run against any service.
