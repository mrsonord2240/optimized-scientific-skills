---
name: ca-bio-causal-genomics-mendelian-randomization
category: Data Analysis
description: Use when estimating exposure-outcome effects from genetic instruments and evaluating MR assumptions.
license: MIT
author: GPTomics
---

# Mendelian randomization

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| Local two-sample summary statistics | `routes/two-sample.md` |
| Outliers or Egger NOME sensitivity | `routes/sensitivity.md` |
| Multivariable MR | `routes/mvmr.md` |
| Overlap, cis, direction or other models | `routes/other.md` |

Assess relevance, independence and exclusion restrictions, ancestry, overlap and pleiotropy. MR is population-level evidence, not individual clinical advice; sensitivity agreement cannot prove assumptions.

Tested source versions include TwoSampleMR 0.6+, MVMR 0.4+, MR-PRESSO 1.0+ and R 4.4+.
