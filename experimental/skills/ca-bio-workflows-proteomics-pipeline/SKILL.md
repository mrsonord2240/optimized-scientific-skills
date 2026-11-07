---
name: ca-bio-workflows-proteomics-pipeline
category: Data Analysis
description: Use when routing search-engine proteomics output to differential abundance.
license: MIT
author: GPTomics
---
# Proteomics workflow
| Intent | Route |
| --- | --- |
| MaxQuant label-free protein matrix | [limma](routes/limma.md) |
| Peptide-level or DIA data | [feature and DIA](routes/feature.md) |
| TMT or SILAC | [labeled](routes/labeled.md) |
| FDR, QC and design caveats | [guidance](references/guidance.md) |
