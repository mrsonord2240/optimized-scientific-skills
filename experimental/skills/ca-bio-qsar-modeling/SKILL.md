---
name: ca-bio-qsar-modeling
description: Use when building or evaluating molecular structure-activity prediction models.
category: Data Analysis
license: MIT
author: GPTomics
---

# QSAR Modeling

Choose a route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Need | Read |
|---|---|
| Chemprop model from assay data | `routes/chemprop.md` |
| Conformal intervals for a sklearn baseline | `routes/conformal.md` |
| Calibration of model probabilities | `routes/calibration.md` |
| Feature attribution or active learning | `routes/analysis.md` |
| Multi-task data | `routes/multitask.md` |
| Related analysis | `routes/related.md` |

Define the endpoint and deployment-relevant split before model selection; fit learned preprocessing on training data only and document applicability and empirical held-out performance.
Source targets: chemprop 2.2.x, scikit-learn >=1.4,<1.6, MAPIE >=0.8,<1.0.
