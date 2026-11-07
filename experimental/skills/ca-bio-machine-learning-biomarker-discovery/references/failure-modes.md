# Biomarker Discovery Failure Modes

Feature-selection failure modes: interpreting minimal-optimal membership, selection-before-CV leakage, the wrong null, winner's curse, single-cell pseudoreplication; reconciliation when selectors disagree.

## Per-Method Failure Modes

### Interpreting minimal-optimal membership as biology
- **Trigger:** Reporting "LASSO selected gene X but not its co-expressed partner Y" as a biological finding.
- **Mechanism:** L1 geometry keeps one vertex of a correlated group arbitrarily; the choice flips across resamples.
- **Symptom:** Selected genes change completely on a different train/test split though accuracy is stable.
- **Fix:** Use elastic net (grouping effect) or report selection *frequencies*; never read membership as importance ordering.

### Selection-before-CV leakage
- **Trigger:** Pick top-k features on all samples, then cross-validate the classifier on those features.
- **Mechanism:** The held-out folds informed which genes were kept; selection is the dominant overfitting capacity in p>>n.
- **Symptom:** Near-perfect CV accuracy, even reproducible on label-permuted (null) data; collapse on an external cohort.
- **Fix:** Selection lives inside the CV fold (Pipeline); estimate by nested CV (machine-learning/model-validation).

### Significance against the wrong null
- **Trigger:** Concluding a signature is real because it significantly predicts outcome.
- **Mechanism:** Random gene sets clear that bar; the transcriptome's proliferation axis is captured by almost any large set (Venet 2011).
- **Symptom:** The signature does not beat a size-matched random signature or a proliferation meta-gene in independent data.
- **Fix:** Benchmark against random-signature and proliferation-meta-gene nulls; require added value over clinical covariates in an *independent* cohort.

### Winner's curse / inflated effect sizes
- **Trigger:** Estimating effect sizes or AUC on the same data used to select features.
- **Mechanism:** Selected features are disproportionately those whose noise inflated their apparent effect (Goring 2001); the inflation can be near-total for small true effects.
- **Symptom:** Discovery AUC much higher than replication; replication is under-powered because it was sized to the inflated effect.
- **Fix:** Estimate effects on an independent split (cross-fitting / data-splitting); size replication for the shrunken effect.

### Pseudoreplication in single-cell selection
- **Trigger:** Treating thousands of cells from a few donors as independent samples.
- **Mechanism:** Cells within a donor are correlated; the effective n is the number of donors.
- **Symptom:** Grossly inflated significance and false discoveries.
- **Fix:** Pseudobulk per donor, select at the donor level (Squair 2021); confront the small true n.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Boruta keeps 200 genes, LASSO keeps 12 | All-relevant vs minimal-optimal answering different questions | Both can be right; pick by goal, do not "average" them |
| A list barely overlaps a published signature | Many disjoint equally-predictive lists exist (Ein-Dor 2005) | Expected, not a contradiction; compare *performance* and stability, not membership |
| High accuracy, low stability index | Resampling accident exploiting a dominant axis | Distrust the specific genes; prefer the lower-accuracy higher-stability candidate |
| FDR-clean list still fails to replicate | FDR controls testing, not selection stability | They are orthogonal; add stability + independent validation |
