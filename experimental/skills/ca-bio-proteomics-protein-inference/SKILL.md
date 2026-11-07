---
name: ca-bio-proteomics-protein-inference
description: Use when grouping protein evidence and estimating protein-level FDR from peptide identifications.
category: Data Analysis
license: MIT
author: GPTomics
---
# Protein inference
| Intent | Route |
| --- | --- |
| Infer groups with OpenMS | [Basic inference](routes/basic.md) |
| Infer Bayesian posteriors | [EPIFANY](routes/epifany.md) |
| Estimate group FDR | [Picked groups](routes/picked.md) |
| Use external CLI methods | [External tools](routes/external.md) |

Retain group membership and decoys through inference; peptide-level FDR does not establish protein-level FDR. See [citations](references/citations.md).
