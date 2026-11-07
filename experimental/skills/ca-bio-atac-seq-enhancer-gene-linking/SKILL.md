---
name: ca-bio-atac-seq-enhancer-gene-linking
category: Data Analysis
description: Use when predicting enhancer-gene links from accessibility and contact data.
license: MIT
author: GPTomics
---

# Enhancer-gene linking

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Data or need | Route |
|---|---|
| Accessibility with optional H3K27ac and contact | `routes/abc.md` |
| ENCODE-rE2G with matched or megamap Hi-C | `routes/re2g.md` |
| Combine thresholded predictions and loop support | `routes/combine.md` |
| Single-cell, co-accessibility or related inputs | `routes/related.md` |

Use calibrated input-specific thresholds, match reference build and document contact source or proxy. Predictions are hypotheses and experimental validation is separate.

Tested source contract: ABC v1.1.2 and main 92ac503; ENCODE-rE2G main d039062.
