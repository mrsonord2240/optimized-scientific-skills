---
name: ca-bio-splicing-quantification
description: Use when measuring splice-site usage or isoform inclusion ratios from short-read RNA sequencing.
category: Data Analysis
license: MIT
author: GPTomics
---

# Splicing Quantification

Choose a route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Data or question | Read |
|---|---|
| rMATS BAM events or existing rMATS results | `routes/rmats.md` |
| Salmon or kallisto transcript TPM | `routes/suppa.md` |
| Annotation-free junction clusters | `routes/leafcutter.md` |
| Complex local splice variations | `routes/majiq.md` |
| Intron retention or microexons | `routes/special.md` |
| Related analyses | `routes/related.md` |

PSI depends on the event and effective-length convention; document inclusion direction, sample order and missing values. For high-confidence reporting reconcile complementary methods.
Source checked: rMATS 4.4.0, SUPPA2 2.4, leafcutter 0.2.9; MAJIQ was not executed.
