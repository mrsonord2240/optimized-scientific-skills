---
name: ca-bio-differential-splicing
description: Use when comparing alternative splicing between conditions using replicated RNA sequencing data.
category: Data Analysis
license: MIT
author: GPTomics
---

# Differential Splicing

Choose a route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Data or question | Read |
|---|---|
| Annotated two-group short-read BAMs | `routes/rmats.md` |
| Novel junctions or annotation-free intron clusters | `routes/leafcutter.md` |
| Existing transcript TPMs, at least four replicates per group | `routes/suppa.md` |
| Few replicates and junction-imbalance correction | `routes/shiba.md` |
| Complex local events or heterogeneous cohorts | `routes/majiq.md` |
| Covariates, multifactor designs or prioritization | `routes/design.md` |
| Related analysis | `routes/related.md` |

Stop for one sample per group: there is no replicate variance. Confirm strandedness, reference and design; perfectly confounded batch and condition cannot be separated. For high-confidence results compare complementary methods, reconcile event definitions and effect directions.
Source checked versions: rMATS 4.4.0, leafcutter 0.2.9, SUPPA2 2.4, Shiba 0.8.2; MAJIQ V3 documentation only.
