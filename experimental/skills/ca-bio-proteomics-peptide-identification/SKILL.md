---
name: ca-bio-proteomics-peptide-identification
description: Use when searching tandem mass spectra and controlling peptide-spectrum match FDR.
category: Data Analysis
license: MIT
author: GPTomics
---
# Peptide identification
| Intent | Route |
| --- | --- |
| Search and rescore | [CLI search](routes/cli.md) |
| Search or filter with pyOpenMS | [pyOpenMS](routes/pyopenms.md) |
| Estimate FDR from tables | [FDR tables](routes/tables.md) |

Record target-decoy strategy, score direction and FDR level. PSM-level FDR does not establish protein-level FDR. See [citations](references/citations.md).
