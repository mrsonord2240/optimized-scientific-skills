Fully confounded batch and condition cannot be separated by a covariate model; redesign or report non-identifiability. Source guidance suggesting otherwise conflicts with its singular-design warning.

## Batch Sources in CRISPR Screens

| Source | Mechanism | Detectable by |
|--------|-----------|---------------|
| Library lot | Different aliquots or PCR amplifications | Gini shift; plasmid-pool sequencing |
| Cell passage cohort | Cells passaged through different periods | PCA Day-0 samples clustering by passage |
| Infection day | Lentivirus titer drifts; FBS lot changes | PCA Day-0 samples cluster by day |
| Cas9 enzyme lot | Cas9 expression heterogeneity | PR-AUC drift across screens |
| Sequencing run | Lane bias, flowcell variant, machine | Per-sample read-count distribution |
| FBS / culture lot | Fetal bovine serum lot variations confound proliferation | Day-0 vs endpoint differential not present in vehicle |
| Tissue-prep batch | In-vivo: animal cohort, surgical day, organ-prep tech | In-vivo screens (see bio-crispr-screens-in-vivo-screens) |

**Critical:** Batch effects in CRISPR screens often correlate with biology (e.g., the drug arm was processed in batch 2 because that's when the drug arrived). This confounds correction. Always check for confounding before applying ComBat.

## Batch Effect Decision Tree

| Diagnostic finding | Recommended correction |
|--------------------|------------------------|
| PCA shows samples cluster by condition, not batch | No correction needed; biology dominates |
| PCA PC1 separates batches, PC2 separates conditions | Apply ComBat with condition passed as `mod` (`references/combat.md`) |
| Batch fully confounded with condition (e.g. all drug in batch 2, all vehicle in batch 1) | Correction will destroy biology; instead redesign next screen with cross-batch balance OR re-analyze with batch in MAGeCK MLE design matrix |
| Day-0 (pre-perturbation) samples cluster by batch | Strong batch effect; ComBat needed (`references/combat.md`) |
| Endpoint samples cluster by batch but not Day-0 | Selection-driven artifact (FBS lot etc); correct or include batch as covariate |
| Replicates within a batch are tight; across-batch much wider | Classic batch effect; ComBat |
| Batch not annotated (unknown technical confounders) | RUV with the NTCs as controls (`references/ruv.md`); or SVA factors as covariates (`references/sva.md`) |
| Each replicate scatters randomly across PCs | Sample-level noise; no batch correction will help |
| Cancer-line panel with multiple batches | Use Chronos (built-in batch and CN modeling). Copy-number bias is a separate, batch-like effect per line: apply CN correction (CRISPRcleanR / Chronos) before batch correction |
| Several screens sharing one library | JACKS (joint efficacy across screens) or Chronos |

## Batch as Explicit Covariate (Preferred for MAGeCK MLE / Chronos)

**Goal:** Model batch and biology in the same regression instead of pre-correcting.

**Approach:** Add batch indicator columns to the MLE design matrix. The fitted beta for condition is the effect after accounting for batch; no pre-correction needed.

```bash
# Design matrix for a screen with 2 batches and 2 conditions
cat > design.txt <<EOF
Samples         baseline    batch2    treatment
Veh_b1_r1       1           0         0
Veh_b1_r2       1           0         0
Drug_b1_r1      1           0         1
Drug_b1_r2      1           0         1
Veh_b2_r1       1           1         0
Veh_b2_r2       1           1         0
Drug_b2_r1      1           1         1
Drug_b2_r2      1           1         1
EOF

mageck mle \
    --count-table counts.txt \
    --design-matrix design.txt \
    --permutation-round 10 \
    --output-prefix batch_aware_mle
```

**Runtime:** at genome scale (~18,000 genes) `mageck mle`'s variance-model permutation is a multi-hour
job by design (a long run is not a hang). For a fast sanity check, run a gene subset or lower
`--permutation-round`; use the full run for the reported result.

**Reproducibility:** the beta estimates are deterministic, but `mageck mle`'s significance comes from
a permutation procedure that is not seeded, so p-values and FDRs move slightly between reruns on
identical input. Fix `--permutation-round` (higher = more stable, linearly slower) and report the
value, or treat borderline FDRs as borderline.

**Why this is preferred:** ComBat shifts counts before testing; the MLE-with-covariates approach correctly propagates uncertainty from the batch term into the condition beta's standard error. ComBat-then-test pretends the corrected counts are noise-free, biasing FDR.

## When NOT to Correct

| Situation | Why correction hurts |
|-----------|----------------------|
| Batch is fully confounded with condition | Correction destroys biology along with batch; redesign or accept |
| Batch effect is smaller than between-replicate noise | Correction adds noise without removing meaningful variance |
| Replicates already correlate >0.95 within and across batches | No batch effect to correct |
| Single-screen analysis | No "batch" to correct; only replicate noise |
| Per-batch sample size <3 | Cannot estimate batch shift reliably; correction is harmful |

## Failure Modes

### ComBat eliminates biological signal

**Trigger:** Batch is correlated with condition (e.g., all drug-arm samples were processed week 2; all vehicle-arm samples week 1).
**Mechanism:** ComBat without a `mod` covariate treats condition variance as batch variance; corrects it away.
**Symptom:** PR-AUC against CEGv2 drops after ComBat correction.
**Fix:** Always supply `mod` (the condition labels); verify by comparing PR-AUC before and after.

### RUV adds noise instead of removing it

**Trigger:** k (number of unwanted factors) set too high.
**Mechanism:** RUV's least-squares decomposition over-fits; "removed" variance includes biology.
**Symptom:** Hits decrease and replicate Pearson drops after correction.
**Fix:** Choose k via cross-validation; default k=1 or 2 for most screens.

### Batch-aware MLE collinear design matrix

**Trigger:** Adding a batch indicator that is fully collinear with another design column (e.g., all of batch 2 is also Day 21).
**Mechanism:** MLE design matrix is singular; betas not estimable.
**Symptom:** MAGeCK MLE errors out or produces NaN betas.
**Fix:** Drop the collinear column; re-design experiment with cross-batch balance.

### ComBat after RUV double-corrects

**Trigger:** Applying multiple corrections sequentially.
**Mechanism:** Both methods remove variance; sequential application removes biology twice.
**Symptom:** All signal gone; counts look uniformly noisy.
**Fix:** Pick one method based on diagnostic; never combine.

### Per-batch sample size too small

**Trigger:** 2 replicates per batch with 3 batches; ComBat estimates batch shift from 2 samples.
**Mechanism:** Insufficient data to estimate batch parameters; high-variance estimates.
**Symptom:** Correction makes some batches worse than uncorrected.
**Fix:** Need ≥3 (preferably 4-6) samples per batch; below this, use covariate modeling instead.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| PC1 batch F vs condition F | F_batch > 10x F_cond -> apply correction | Standard variance-decomposition diagnostic |
| ComBat min samples per batch | ≥3, ideally 4-6 | Empirical Bayes prior estimation |
| RUV `k` (unwanted factors) | k=1 default; k=2 if multiple known batch sources | Risso 2014; cross-validate |
| NTCs needed for NTC-anchored norm | ≥500 in library | Stable median |
| Post-correction PCA check | Batches must overlap in PC1/PC2 plot | Visual sanity check |
| Post-correction PR-AUC | Should be same or higher than pre | If lower, correction destroyed biology |

## Validation Checklist

After applying correction:

- [ ] Corrected matrix has no NaN/Inf values (ComBat can return all-NaN with exit code 0)
- [ ] Features in `uncorrected` are flagged or excluded in hit calling
- [ ] PCA: batches now overlap (visual)
- [ ] Within-batch Pearson preserved (should be unchanged)
- [ ] Across-batch Pearson improved
- [ ] CEGv2 PR-AUC preserved or higher
- [ ] NTC distribution stable across batches
- [ ] No new outlier samples introduced
- [ ] Hit list compared with the uncorrected hit list; every difference explained by the batch effect, not lost biology

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| PR-AUC drops after ComBat | Batch confounded with condition | Add `mod` covariate; or redesign |
| `combat_correct()` reports N features uncorrected | No variance left after batch and condition (e.g. all-zero guides), so ComBat cannot fit them | Expected; flag or exclude the returned `uncorrected` index in hit calling, or model batch as a covariate |
| MAGeCK MLE NaN beta after adding batch column | Collinear design matrix | Drop collinear column |
| Replicates still cluster by batch after RUV | k too low | Increase k; cross-validate |
| Replicates lose internal cohesion after correction | Over-correction | Reduce k or revert |
| NTC-anchored norm worse than median | Too few NTCs | Use median; add NTCs to next library |
| Sequencing-run-level batch survives ComBat | Non-linear sequencing effect | Pre-normalize with `mageck count --norm-method control` first |
