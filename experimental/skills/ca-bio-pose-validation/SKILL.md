---
name: ca-bio-pose-validation
description: Use when assessing physical plausibility and strain of docked or generated protein-ligand poses.
category: Data Analysis
license: MIT
author: GPTomics
---

# Pose Validation

Choose a route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Need | Read |
|---|---|
| PoseBusters checks on an SDF and receptor | `routes/posebusters.md` |
| Relative ligand strain | `routes/strain.md` |
| Batch pose filtering | `routes/batch.md` |
| Aromatic planarity diagnostic | `routes/planarity.md` |
| Related analysis | `routes/related.md` |

Inspect actual checks and loading statuses; configuration and version determine which tests run. A plausibility pass does not establish a correct binding mode.
Source checked: PoseBusters 0.6+, RDKit 2024.09+.
