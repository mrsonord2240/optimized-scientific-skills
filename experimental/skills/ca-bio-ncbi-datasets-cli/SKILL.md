---
name: ca-bio-ncbi-datasets-cli
description: Use when retrieving genome, gene, ortholog, or virus data with NCBI Datasets CLI.
license: MIT
author: GPTomics
category: Data Analysis
---
# NCBI Datasets CLI

| Intent | Route |
|---|---|
| Genome download, summary, or Python wrapper | [Genome](routes/genome.md) |
| Bulk discovery and transfer | [Bulk](routes/bulk.md) |
| Gene and ortholog metadata | [Gene](routes/gene.md) |
| Virus genome and metadata | [Virus](routes/virus.md) |

This experimental normalization preserves source examples checked against CLI 18.37.0; verify the installed CLI and field names before use. See [background](references/background.md) and [citations](references/citations.md). Raw SRA reads are outside this skill; use `bio-sra-data`. For tree-aware orthology use `bio-ortholog-inference`.
