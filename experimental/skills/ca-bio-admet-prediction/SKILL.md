---
name: ca-bio-admet-prediction
description: Use when estimating small-molecule ADMET properties or drug-likeness for lead triage.
category: Data Analysis
license: MIT
author: GPTomics
---

# ADMET Prediction

Choose a route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Need | Read |
|---|---|
| Offline ADMET-AI batch with input gate | `routes/offline.md` |
| Imported hosted ADMETlab results | `routes/hosted.md` |
| Custom endpoint model and calibration | `routes/custom.md` |
| Physicochemical properties or structural alerts | `routes/properties.md` |
| hERG, CYP or model disagreements | `routes/endpoints.md` |
| Related analysis | `routes/related.md` |

Predictions support prioritization, not clinical or regulatory decisions. Evaluate the applicability domain, preserve flagged inputs and validate important results experimentally.
Source checked: ADMET-AI 2.0.1, chemprop 2.3.1, RDKit 2024.09+.
