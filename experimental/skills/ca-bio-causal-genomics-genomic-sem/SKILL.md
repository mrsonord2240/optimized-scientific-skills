---
name: ca-bio-causal-genomics-genomic-sem
category: Data Analysis
description: Use when modeling shared genetic factors or multivariate GWAS effects from summary statistics.
license: MIT
author: GPTomics
---

# Genomic SEM

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| LDSC covariance and factor model | `routes/factor.md` |
| Common-factor GWAS and heterogeneity | `routes/gwas.md` |
| Custom or exploratory models | `routes/custom.md` |
| Stratified factor heritability | `routes/stratified.md` |
| Compare MTAG | `routes/mtag.md` |

Use the full LDSC S and sampling covariance V, match reference ancestry, and check identification and model fit. Results concern population-level architecture, not individual risk or diagnosis.

Tested source pairing: GenomicSEM 0.0.5 with lavaan 0.6.19; source reports incompatibility with lavaan >= 0.7.
