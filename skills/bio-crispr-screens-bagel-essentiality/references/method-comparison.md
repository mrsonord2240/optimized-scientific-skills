# BAGEL2 vs MAGeCK vs drugZ

> Moved out of SKILL.md. "Below" and "Failure Modes" in the text refer to sections of SKILL.md.

## Comparing BAGEL2, MAGeCK, drugZ

| Property | BAGEL2 | MAGeCK | drugZ |
|----------|--------|--------|-------|
| Statistical framework | Bayes factor with reference sets | NB GLM | Bidirectional Z-score |
| Calibrated against | CEGv2 / NEGv1 | Internal null | Vehicle distribution |
| Tumor suppressor detection | YES | Limited (RRA positive-selection score) | YES |
| Best for | Essentiality classification | General hit calling | Chemogenomic drug screens |
| Output | Bayes factor + CI | FDR + LFC | Z-score + FDR per direction |
| Hit threshold | BF >6 | FDR <0.05 | FDR <0.05 |
| Library calibration | Indirect (reference set) | None | None |

**Pick by case:** drug-modifier screen wanting both directions -> BAGEL2 or drugZ; multi-cell-line panel -> Chronos (DepMap) rather than BAGEL2; custom library with <4 sgRNAs/gene or a non-cancer line -> see Failure Modes.

**Reconciliation:** BF >6 ≈ 90% posterior probability (Hart 2017 G3); commonly treated as roughly MAGeCK FDR 0.05 by convention. BAGEL2 hits absent from MAGeCK suggest weak signal that BAGEL2's reference anchoring detects but MAGeCK's null-based test misses; verify by inspecting per-sgRNA contributions.
