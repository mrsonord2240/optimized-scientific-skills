# Design guidance

A missing factor-by-batch cell is not automatically perfect aliasing; the assignment helper uses a conservative complete-coverage rule. Assess full design matrix rank.

## The Single Most Important Modern Insight -- No Post-Hoc Method Recovers a Confounded Design

When a technical factor is perfectly aliased with the biological factor (all treated in batch 1, all controls in batch 2), batch and condition occupy the same column space and are mathematically non-identifiable; ComBat, SVA, and RUV cannot separate them, and "removing the batch effect" removes the biology with it. Worse, on *partially* confounded or merely **unbalanced** designs, mean-centering batches can *manufacture* false positives and inflate downstream confidence — Nygaard, Rødland & Hovig 2016 *Biostatistics* 17:29 showed a pipeline that returned >1000 spurious DE probes where the honest analysis (batch kept *in the model*) found 11. The operative rules: (1) balance the biological variable across batches **at design time** because correction is not a rescue (Leek 2010 *Nat Rev Genet* 11:733); (2) for inference, keep batch *in the model* so its degrees of freedom are charged honestly — reserve a batch-"cleaned" matrix for visualization and clustering only. The experimental unit, not the measurement, still defines replication (see experimental-design/randomization-blocking).

## Confounding vs Blocking vs Nuisance -- the Causal-Graph View

Batch is the genomics face of confounding, and not all metadata should be "adjusted for". A **confounder** is a common cause of treatment and outcome (adjust for it); a **mediator** lies on the causal path (adjusting removes signal); a **collider** is a common effect (adjusting *induces* spurious association). "Adjust for everything measured" is therefore wrong in general — conditioning on a collider opens a backdoor path. In a properly *randomized* design the biological variable has no confounders by construction, so batch is handled by **balanced assignment + a block/covariate term**, not by scrubbing every measured variable.

## Algorithmic Taxonomy -- Design Strategies for Technical Variation

| Strategy | What it does | When to use | Fails when |
|----------|--------------|-------------|------------|
| Balanced (orthogonal) assignment | every condition appears equally in every batch | always achievable when batches hold >=1 of each condition | not all conditions fit per batch |
| Block randomization across batches | randomize condition within each batch | batch = a block; conditions fit per batch | batch variance is genuinely zero (rare) |
| Incomplete block + batch in model | conditions split across smaller batches, batch term retained | plate/chip smaller than #conditions | unbalanced split inflates artifacts (Nygaard 2016) |
| Reference / bridge sample per batch | shared anchor measured in every batch | cross-batch normalization (TMT proteomics, large cohorts) | anchor not representative |
| Multiplexing + demultiplexing | pool biological units in one lane, split by barcode/genotype | breaking the donor<->lane confound (scRNA-seq) | insufficient SNPs/hashes to assign |
| Run-order randomization | randomize processing/injection order | position/time gradients (LC-MS, plate edge) | order set by convenience |

## Decision Tree by Scenario

| Scenario | Recommended design | Why |
|----------|--------------------|-----|
| 24 samples, 3 batches, 2 conditions | balanced: 4 of each condition per batch | batch orthogonal to condition; estimable |
| Conditions outnumber batch capacity | incomplete block; keep batch in the DE model | preserves estimability; no scrubbing |
| Large cohort across many runs | include a shared reference sample per batch | enables cross-batch normalization |
| scRNA-seq, several donors, few lanes | pool donors per lane, demultiplex (demuxlet / hashing) | removes donor<->lane confound (Kang 2018) |
| Hidden/unknown technical structure suspected | estimate surrogate variables (SVA), include in model | captures unmodeled variation (Leek & Storey 2007) |
| Design already confounds batch with condition | redesign; no correction will rescue it | non-identifiable (Nygaard 2016; Leek 2010) |
| General randomization / blocking / unit choice | -> experimental-design/randomization-blocking | foundational design structure |
| Running ComBat-seq / RUVSeq / SVA on real data | -> differential-expression/batch-correction | execution lives there; this skill decides |

## Downstream Correction -- Choose by Design, Execute Elsewhere

Correction method selection is a design decision; the execution lives in differential-expression/batch-correction (and single-cell/batch-integration for scRNA-seq). Prefer keeping batch in the analysis model over producing a "cleaned" matrix for inference.

| Method | When it applies | Assumption / caveat | Owner of execution |
|--------|-----------------|---------------------|--------------------|
| Batch as a model covariate | batch known, balanced | charges df honestly; the default for inference | differential-expression |
| ComBat-seq | known batch, integer counts | batch ~ orthogonal to biology; Nygaard caveat if unbalanced | differential-expression/batch-correction |
| ComBat (parametric eB) | known batch, log-normalized data | Gaussian; not for raw counts | differential-expression/batch-correction |
| RUVSeq (RUVg/RUVs/RUVr) | negative-control genes/samples available | controls must be truly null to the biology | differential-expression/batch-correction |
| SVA | hidden/unknown structure | surrogate variables can absorb biology if confounded | this skill estimates; DE consumes |
| limma removeBatchEffect | visualization/clustering ONLY | not for the hypothesis test | data-visualization |
| Harmony / scVI / Seurat anchors | scRNA-seq integration | integration, not DE inference | single-cell/batch-integration |

## Detecting Hidden Batch Effects (SVA)

**Goal:** Estimate unmodeled technical structure (hidden batches) so it can be included in the downstream model.

**Approach:** Fit a model matrix for the biological variable and a null matrix, estimate the number of surrogate variables, then compute them for inclusion as covariates in the DE analysis. `sva()`/`num.sv()` require a **complete matrix** -- they stop with `infinite or missing values in 'x'` on any NA/NaN/Inf, which is the normal state of an LC-MS or LFQ intensity matrix before imputation. Decide explicitly how to handle that rather than letting the error surface unexplained or silently dropping rows without saying so, and report how many cells were non-finite:

- **Option A** -- keep only complete-observation features. Simple and defensible for SV estimation, but it biases the surrogate variables toward abundant, well-detected features; confirm that is acceptable before trusting the SVs for low-abundance biology.
- **Option B** -- impute first (normalization-qc covers QRILC/missForest for left-censored proteomics/metabolomics intensities), then run `sva()` on the imputed matrix. The imputer shapes the estimated SVs, so record which one was used.

`run_sva_safely()` in `examples/batch_design.R` implements Option A end to end (count and report the non-finite cells, restrict to complete features, assert finiteness, then `num.sv()` + `sva()`); copy it. The core calls, once the matrix is complete:

```r
library(sva)
mod  <- model.matrix(~ condition, data = colData)   # full model
mod0 <- model.matrix(~ 1, data = colData)           # null model
n_sv <- num.sv(expr_normalized, mod)                # estimate number of hidden batches (complete matrix only)
svobj <- sva(expr_normalized, mod, mod0, n.sv = n_sv)
# Add svobj$sv to the design used by differential-expression/de-results; do NOT subtract them
# from the data for the test (subtracting is for visualization only).
```

## Per-Method Failure Modes

### Batch confounded with condition
- **Trigger:** all of one condition processed in one batch/run.
- **Mechanism:** batch and condition are aliased -> non-identifiable.
- **Symptom:** condition effect vanishes (or an artifact appears) after correction.
- **Fix:** redesign with balanced assignment; no post-hoc method recovers it (Leek 2010).

### ComBat on unbalanced groups
- **Trigger:** ComBat applied when condition is partially confounded with batch.
- **Mechanism:** mean-centering batches injects between-group differences and understates residual variance.
- **Symptom:** inflated DE counts and over-confident downstream inference (Nygaard 2016).
- **Fix:** keep batch in the model (ComBat-seq with the biological covariate, or batch as a DE covariate); best is balanced design.

### Subtracting surrogate variables before testing
- **Trigger:** feeding an SV-"cleaned" matrix into the DE test.
- **Mechanism:** double-counts the adjustment and loses degrees of freedom.
- **Symptom:** anti-conservative p-values.
- **Fix:** include SVs as covariates in the model; reserve cleaned matrices for plots.

### Adjusting for a collider
- **Trigger:** "adjust for everything measured" includes a downstream/common-effect variable.
- **Mechanism:** conditioning on a collider opens a spurious path.
- **Symptom:** associations that appear only after adjustment.
- **Fix:** adjust for confounders (common causes), not mediators or colliders.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Balance every condition equally across batches | Leek 2010 *Nat Rev Genet* 11:733 | makes batch estimable and removable |
| Unbalanced ComBat can inflate DE (>1000 vs 11 in one case) | Nygaard 2016 *Biostatistics* 17:29 | mean-centering injects group differences |
| ~50 SNPs/cell suffice to demultiplex pooled donors | Kang 2018 *Nat Biotechnol* 36:89 | breaks donor<->lane confound |
| Keep batch in the model for inference; clean only for viz | Nygaard 2016; Leek 2010 | honest degrees of freedom |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `sva()`/`num.sv()` error: infinite or missing values in 'x' | proteomics/metabolomics matrix has NAs (normal before imputation) | filter to complete features or impute first, then re-run (see Detecting Hidden Batch Effects above) |
| Condition effect disappears after ComBat | batch confounded with condition | balance at design time |
| Inflated DE list after batch correction | unbalanced ComBat | keep batch in model; ComBat-seq with covariate |
| scRNA-seq donor effect equals lane effect | one donor per lane | pool + demultiplex (demuxlet / hashing) |
| Spurious associations after "adjusting for all metadata" | conditioning on a collider/mediator | adjust only for confounders |
| Cannot reconstruct who/when/which-lot | no metadata captured | record date, lot, operator, lane, position |

## Reproducibility Metadata

Record for every sample, because these become the batch/blocking variables: processing date, reagent and kit lot numbers, operator, instrument/flow-cell/lane and well/plate position, library prep batch, and any protocol deviations. Unrecorded technical variation cannot be modeled or balanced after the fact, and is a leakage source for any downstream machine-learning model (see machine-learning/model-validation).

## References

- Leek JT, Scharpf RB, Bravo HC, Simcha D, Langmead B, Johnson WE, Geman D, Baggerly K, Irizarry RA. 2010. Tackling the widespread and critical impact of batch effects in high-throughput data. *Nat Rev Genet* 11:733-739.
- Nygaard V, Rødland EA, Hovig E. 2016. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. *Biostatistics* 17:29-39.
- Johnson WE, Li C, Rabinovic A. 2007. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics* 8:118-127.
- Zhang Y, Parmigiani G, Johnson WE. 2020. ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genom Bioinform* 2:lqaa078.
- Leek JT, Storey JD. 2007. Capturing heterogeneity in gene expression studies by surrogate variable analysis. *PLoS Genet* 3:e161.
- Gagnon-Bartsch JA, Speed TP. 2012. Using control genes to correct for unwanted variation in microarray data. *Biostatistics* 13:539-552.
- Kang HM, Subramaniam M, Targ S, et al. 2018. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. *Nat Biotechnol* 36:89-94.
