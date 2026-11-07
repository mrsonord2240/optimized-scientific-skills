---
name: ca-bio-causal-genomics-fine-mapping
category: Data Analysis
description: Use when fine-mapping association signals into candidate causal variants and credible sets.
license: MIT
author: GPTomics
---

# Fine-mapping

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| Summary statistics and signed LD | `routes/susie.md` |
| Independent FINEMAP or DAP-G analysis | `routes/independent.md` |
| Multiple ancestries | `routes/susiex.md` |
| Functional priors | `routes/priors.md` |
| Colocalization, HLA, or other models | `routes/related.md` |

Harmonize alleles and SNP order with ancestry-matched LD. Report LD diagnostics, credible-set size and purity, coverage, PIPs, sample size and prior. Credible variants are candidates, not proven causal.

Tested source versions: susieR 0.14.2, coloc 5.2.3, FINEMAP 1.4.2, SuSiEx 1.1.2.
