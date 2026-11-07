---
name: ca-bio-causal-genomics-colocalization-analysis
category: Data Analysis
description: Use when testing whether association signals share causal variants across traits or molecular QTLs.
license: MIT
author: GPTomics
---

# Colocalization

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| Single-causal GWAS and QTL | `routes/abf.md` |
| Multiple signals with signed LD | `routes/susie.md` |
| Conditional signal pairs or per-SNP CLPP | `routes/multicausal.md` |
| Many traits or multiple omics | `routes/multitrait.md` |
| SMR and HEIDI | `routes/smr.md` |
| MHC, inversion, or unsupported method | `routes/exception.md` |

Harmonize alleles and genome build before analysis; assess power and ancestry. Report all posterior hypotheses and prior sensitivity, with PP.H3 alongside PP.H4. Do not interpret standard single-causal results in the MHC or chr8 inversion.

Tested source versions: coloc 5.2.3+, susieR 0.12.35+, R >= 4.1.
