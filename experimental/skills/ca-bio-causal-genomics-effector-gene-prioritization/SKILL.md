---
name: ca-bio-causal-genomics-effector-gene-prioritization
category: Data Analysis
description: Use when prioritizing candidate effector genes at GWAS loci from multiple evidence streams.
license: MIT
author: GPTomics
---

# Effector-gene prioritization

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| Integrate per-gene evidence | `routes/concordance.md` |
| MAGMA gene association | `routes/magma.md` |
| PoPS similarity ranking | `routes/pops.md` |
| Open Targets L2G | `routes/l2g.md` |
| cS2G lookup | `routes/cs2g.md` |
| FLAMES pretrained scorer | `routes/flames.md` |
| Fine-mapping, coloc and chromatin evidence | `routes/related.md` |

These are population-level research hypotheses, not individual clinical determinations. Match tissue, genome build and ancestry; allow multiple candidate genes and distinguish missing from negative evidence.

Tested source versions include MAGMA 1.10, R 4.3 and Python 3.9-3.11; external API schemas may change.
