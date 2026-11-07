---
name: ca-bio-causal-genomics-mediation-analysis
category: Data Analysis
description: Use when estimating direct and indirect effects through molecular or other mediators.
license: MIT
author: GPTomics
---

# Mediation analysis

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Need | Route |
|---|---|
| Observational single mediator | `routes/observational.md` |
| Exposure-mediator interaction | `routes/four-way.md` |
| High-dimensional mediator screening | `routes/hima.md` |
| Double-ML mediation | `routes/dml.md` |
| MR or longitudinal mediation | `routes/other.md` |

Estimates rely on a justified causal model and untestable assumptions; assess sensitivity to unmeasured confounding. Results are population-level research evidence, not individual treatment guidance.

Tested source versions include mediation 4.5.0+, HIMA 2.3.4, CMAverse 0.1.0+ and causalweight 1.1.4.
