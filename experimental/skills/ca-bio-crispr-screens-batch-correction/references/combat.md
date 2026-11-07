## ComBat Empirical-Bayes Correction

**Goal:** Remove batch-specific location and scale shifts while preserving biological condition signal.

**Approach:** Log-transform counts, fit ComBat with the condition passed as `mod` (so the model knows which signal to preserve), back-transform. Features ComBat cannot fit are left as raw counts and returned by name.

```bash
python scripts/combat_correct.py counts.txt metadata.txt corrected.txt uncorrected.txt \
    --batch-col batch --condition-col condition     # --condition-col "" runs without a covariate
```

`counts.txt` is the tab-separated count table (first column sgRNA id; Gene and other non-sample columns pass through), `metadata.txt` has the sample name in column one plus the batch and condition columns. From Python: `from combat_correct import combat_correct` (with `scripts/` on the path) returns `(corrected_counts, uncorrected)`. The script log2-transforms, drops the features ComBat cannot fit, runs `pycombat` with the condition as `mod`, back-transforms, and returns the dropped features as raw counts. Its docstring explains why: a feature with no variance left after batch and condition are regressed out corrupts ComBat's shared prior (all-NaN output or a silently constant feature). The pre-fit filter prevents this, the NaN check after the fit is only a backstop, so do not remove the filter.

Run `combat_correct()` only when the diagnostic says batch dominates; on a batch-free design ComBat has nothing to remove and adds noise. Keep the `uncorrected` index and flag those features (or exclude them) in hit calling: they still carry their batch effect. Model them with batch as a covariate instead ("Batch as Explicit Covariate" in `SKILL.md`).

**Critical caveat:** ComBat assumes batch effects are linear shifts of mean and variance in log space. Non-linear effects (e.g., gene-specific batch sensitivity) remain. Always re-check PCA after correction to confirm batches now overlap.
