---
name: ca-bio-atac-seq-deep-learning-atac
category: Data Analysis
description: Use when modeling ATAC-seq accessibility from sequence or interpreting trained sequence models.
license: MIT
author: GPTomics
---

# ATAC sequence models

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| Train bias-corrected accessibility model | `routes/train.md` |
| Score variants using a trained model | `routes/variants.md` |
| Motif discovery from attributions | `routes/motifs.md` |
| Long-range track variant effects | `routes/enformer.md` |
| Per-cell accessibility embeddings | `routes/scbasset.md` |
| Classical ATAC analysis or related tasks | `routes/related.md` |

Use separate compatible environments as described in `references/environments.md`. Apply the model QC gate before variant or motif interpretation; match genome and cell type.

Tested contract: chromBPNet 1.0.1, bpnet-lite 1.0.0, tangermeme 1.5.0, modisco-lite 2.4.0 and enformer-pytorch 0.8.12.
