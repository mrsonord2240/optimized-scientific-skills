# Reading the output

What a z-score is, per motif and sample: raw deviation is (observed - expected) / expected over motif-containing peaks; the bias-corrected deviation subtracts the mean of matched-background peak sets (matched on GC and mean accessibility); the z-score divides by the background SD. Positive means more accessible than the population average. `deviations()` returns bias-corrected deviations, `deviationScores()` the z-scores.

| Variability | Typical z-score range | Reading |
|-------------|-----------------------|---------|
| under 1 | -1 to +1 | Roughly constant |
| 1-2 | -2 to +2 | Modest; condition-driven possible |
| 2-5 | -3 to +5 | Strong cross-sample variability |
| over 5 | -5 to +10 | Major driver, or a sparse-data artefact (check peak count and depth first) |

The ranges are inherited guidance, not measured, and depend on the contrast.

- chromVAR asks whether peaks containing a motif are more or less accessible than expected, and works when the motif has well over 100 sites. Footprinting asks whether one site is bound (`routes/other-skills.md`). They complement each other.
- Top motifs that disagree between tools usually reflect different motif databases; z-scores that correlate but differ in size reflect different background sampling. Do not compare stored z-scores from different sources.
- Decoupling chromVAR z-scores into TF activity (for example with decoupleR) is not valid: motif z-scores are not gene targets. TF-to-target inference needs paired RNA (SCENIC+).
- Sequence-model alternatives (scBasset, Enformer-derived activity) are not part of this Skill; see atac-seq/deep-learning-atac.

Done when the answer states the contrast, the seed and database, and whether sparsity could explain the top motifs.
