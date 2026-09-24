---
name: bio-experimental-design-power-analysis
description: Calculates statistical power for high-dimensional genomics experiments (bulk RNA-seq, scRNA-seq, ATAC-seq, ChIP-seq, methylation, proteomics) under negative-binomial count models using RNASeqPower, PROPER, and simulation via powsimR, distinguishing per-gene from marginal (transcriptome-wide) power, the role of mean expression and dispersion, and the sequencing-depth-versus-replicate tradeoff. Covers simulation as the honest default for overdispersed counts, FDR-aware average power versus single-test power, observed/post-hoc power as an anti-pattern, and the winner's-curse / Type-S / Type-M consequences of underpowering. Use when planning replicate number for a sequencing experiment, deciding whether to add depth or samples, choosing closed-form versus simulation power, estimating power from pilot dispersions, or justifying replication in a grant. For clinical-trial power see clinical-biostatistics/power-and-sample-size; for the inverse sample-size question see experimental-design/sample-size.
tool_type: r
primary_tool: RNASeqPower
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: RNASeqPower 1.46.0, PROPER 1.38.0, DESeq2 1.46.0, edgeR 4.4.2, pwr 1.3.0 (checked 2026-09-17). powsimR 1.2+ (GitHub) is named below as an optional, heavier alternative for scRNA-seq; it is not required for any code block in this Skill.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws an error, introspect the installed package and adapt to the actual API. Notes: `RNASeqPower::rnapower()` solves for whichever of `n` or `power` is omitted, and despite its name applies to any negative-binomial per-feature count assay, not only RNA-seq (its own vignette is titled "Sample Size for RNA-Seq and **similar** Studies") — see the ATAC/ChIP/methylation section below; PROPER is a multi-step pipeline (`RNAseq.SimOptions.2grp` -> `simRNAseq` -> `runSims` -> `comparePower`) whose `RNAseq.SimOptions.2grp` hard-codes `sim.seed = 11111` when `sim.seed` is not supplied, which is why repeated runs of the blocks below are byte-identical without an explicit `set.seed()`; pass `sim.seed = <n>` explicitly for a documented, intentional re-run. powsimR, if installed separately, is GitHub-only and its `estimateParam`/`Setup`/`simulateDE` signatures drift — pin a commit SHA for reproducible work. Verify each against the installed help before relying on argument names.

**Parameter ranges (checked on RNASeqPower 1.46.0, pwr 1.3.0):** `alpha` and `power` must lie in (0, 1); `effect` is a positive fold change other than 1 (`0.5` gives the same result as `2`, so depletion needs no separate handling); `cv`, `depth` and `n` must be positive. `pwr.t.test` stops with an error on `power`/`sig.level` outside [0, 1], but `rnapower()` does not validate: `power = 1.2` returns `NaN` (with a warning), `power = 1` or `effect = 1` returns `Inf`, `alpha = 1.5` returns power 0.999 and `effect = 0` returns power 1 with no warning, and a negative `cv` is silently treated as positive. Check inputs before calling, and treat any `NaN`/`Inf`/implausible result as a bad input rather than a finding.

# Power Analysis for Genomics Experiments

**"How many replicates does my sequencing experiment need?"** -> Compute the probability of detecting a biologically meaningful effect given replicate number, sequencing depth, and biological variability — modeling counts as negative-binomial and recognizing that power is a per-gene quantity, not one number for the whole transcriptome.
- R: `RNASeqPower::rnapower()` — closed-form NB power/sample size; `PROPER`, `powsimR` — simulation from the mean-dispersion trend

## The Single Most Important Modern Insight -- Genomics Power Is Per-Gene; Simulate, and Never Report Observed Power

Power in a sequencing experiment is not a single number. It is a per-gene quantity that depends on that gene's mean expression and dispersion, so the honest summary is the **marginal (average) power** across the expression distribution at a target FDR — the expected discovery rate. A single coefficient of variation plugged into a closed-form formula mis-states power for low- and high-expressed genes alike, because dispersion varies systematically with the mean; the defensible default for count data is **simulation from the empirical mean-dispersion trend** (PROPER, Wu 2015 *Bioinformatics* 31:233; powsimR, Vieth 2017 *Bioinformatics* 33:3486). The second rule is negative: **observed (post-hoc) power is information-free.** Computed from the effect a study actually estimated, it is a one-to-one function of the p-value and cannot explain a null result (Hoenig & Heisey 2001 *Am Stat* 55:19). Power is a design-stage quantity, computed for hypothesized effects before data exist. Underpowering does not merely miss true effects — it makes the significant ones overstate magnitude (Type-M) and sometimes reverse sign (Type-S), lowering the chance a significant call is real (Button 2013 *Nat Rev Neurosci* 14:365; Gelman & Carlin 2014 *Perspect Psychol Sci* 9:641).

## Algorithmic Taxonomy

| Approach | Model | Tool | Strength | Fails / costs when |
|----------|-------|------|----------|--------------------|
| NB closed-form | negative-binomial, single CV/dispersion | `RNASeqPower::rnapower` | fast; transparent; grant-ready | one CV cannot represent the mean-dispersion trend |
| Simulation, parametric | NB with mean-dispersion relationship | `PROPER` | honest marginal power + EDR at target FDR | needs a dispersion model / pilot |
| Pseudobulk NB | donor-level counts, cells aggregated first | `edgeR` (or `RNASeqPower`/`PROPER` on the aggregated matrix) | scRNA-seq population DE; no extra dependency | discards within-donor structure; needs enough donors |
| Simulation, empirical | resampled from pilot (incl. dropout) | `powsimR` (optional, GitHub-only) | bulk AND scRNA-seq; realistic; models dropout directly | GitHub-only; heavier; version drift; not required — pseudobulk NB covers the population-power question |
| Gaussian closed-form | t-test / Cohen's d | `pwr::pwr.t.test` with proteome-wide multiplicity correction | per-feature ATAC/proteomics after transform | wrong for raw counts; ignores overdispersion; understates n by 3-7x if multiplicity is skipped |
| Effect-inflation design analysis | retrodesign for Type-S/Type-M | `retrodesign` (Gelman) | exposes exaggeration in noisy small-n | needs a plausible true effect |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Bulk RNA-seq, pilot data available | PROPER/powsimR simulation from pilot dispersions | matches the real mean-dispersion trend |
| Bulk RNA-seq, no pilot, quick grant number | `rnapower()` with a literature CV, stated as approximate | transparent; flag as conservative-to-rough |
| scRNA-seq cross-condition DE | aggregate to pseudobulk (edgeR), power on donor count — see scRNA-seq section | population power is set by donors, not cells |
| ATAC/ChIP/methylation per-region | `rnapower()`/PROPER directly on per-region counts — see ATAC/ChIP/methylation section | same NB machinery as RNA-seq; overdispersed counts; per-region power |
| Proteomics (continuous, log-abundance) | `pwr::pwr.t.test` per protein with proteome-wide multiplicity correction — see Proteomics section | Gaussian after transform; MNAR matters; raw per-protein alpha understates n 3-7x |
| Justifying a null result post-hoc | report CI / effect size, NOT observed power | post-hoc power is uninformative (Hoenig-Heisey) |
| Fixed budget: depth vs replicates | favor replicates past ~10-20M mapped reads | biological variance dominates (Liu 2014) |
| Clinical-trial endpoint | -> clinical-biostatistics/power-and-sample-size | regulated regime, different machinery |

## Closed-Form NB Power -- RNASeqPower

**Goal:** Get a fast, transparent power or replicate number for bulk RNA-seq from depth, biological CV, and fold change.

**Approach:** Supply per-gene depth, biological coefficient of variation, the fold change to detect, and alpha; supply `n` to get power, or `power` to get the required `n`. Use the checked wrapper below rather than calling `rnapower()` directly: RNASeqPower can return plausible numbers for invalid values. Treat the result as a single-gene approximation and sanity-check against simulation.

```r
library(RNASeqPower)

# Exactly one of n, power, or effect may be omitted for RNASeqPower to solve.
# This wrapper intentionally accepts scalar prospective-design inputs only.
validate_rnapower_inputs <- function(depth, cv, effect = NULL, alpha = 0.05,
                                     n = NULL, power = NULL) {
  check_scalar <- function(x, name, predicate, rule) {
    if (length(x) != 1L || !is.finite(x) || !predicate(x)) {
      stop(sprintf("%s must be one finite scalar %s", name, rule), call. = FALSE)
    }
  }
  check_scalar(depth, "depth", function(x) x > 0, "> 0")
  check_scalar(cv, "cv", function(x) x > 0, "> 0")
  check_scalar(alpha, "alpha", function(x) x > 0 && x < 1, "in (0, 1)")
  if (!is.null(n)) check_scalar(n, "n", function(x) x > 0, "> 0")
  if (!is.null(power)) check_scalar(power, "power", function(x) x > 0 && x < 1, "in (0, 1)")
  if (!is.null(effect)) check_scalar(effect, "effect", function(x) x > 0 && x != 1, "> 0 and != 1")
  if (sum(vapply(list(n, power, effect), is.null, logical(1))) != 1L) {
    stop("Provide exactly two of n, power, and effect; RNASeqPower solves the omitted quantity.", call. = FALSE)
  }
  invisible(TRUE)
}

checked_rnapower <- function(depth, cv, effect = NULL, alpha = 0.05,
                             n = NULL, power = NULL) {
  validate_rnapower_inputs(depth, cv, effect, alpha, n, power)
  args <- Filter(Negate(is.null), list(depth = depth, n = n, cv = cv,
                                        effect = effect, alpha = alpha, power = power))
  do.call(RNASeqPower::rnapower, args)
}

# depth = per-gene coverage (see "Depth Units" below, NOT total library size);
# cv = biological coefficient of variation; effect = fold change
checked_rnapower(depth = 20, n = 5, cv = 0.4, effect = 2, alpha = 0.05)          # solves for POWER
checked_rnapower(depth = 20, cv = 0.4, effect = 2, alpha = 0.05, power = 0.80)   # solves for n per group
```

## Depth Units -- Converting a Real Read Budget to `depth`

**Goal:** Turn a stated sequencing budget (e.g. "20 million reads per sample") into the `depth` argument `rnapower()` expects, instead of guessing.

`depth` is **not** total library size — it is the average per-gene coverage, and RNASeqPower's own vignette (Hart et al. 2013, `samplesize.Rnw`, installed with the package: `vignette('samplesize', package = 'RNASeqPower')`) gives the conversion: across the studies they examined, **85-95% of targets had coverage >= 0.1 per million mapped reads** — i.e. a 40-million-read library gives `depth ~= 4` (`40 * 0.1`) for the majority of genes. This is a conservative floor most genes clear, not their mean coverage; well-expressed genes run far higher. Use it to size `depth` from a real budget rather than reusing the vignette's illustrative `depth = 20` (which implies ~200M mapped reads/sample and is far deeper than a typical bulk RNA-seq budget):

```r
library(RNASeqPower)
reads_millions <- 20                       # e.g. 20M mapped reads/sample
depth_conservative <- 0.1 * reads_millions # majority-of-targets floor, per the package vignette
cat('20M reads/sample -> depth ~=', depth_conservative, '(conservative, per Hart 2013)\n')
checked_rnapower(depth = depth_conservative, n = 14, cv = 0.3, effect = 1.5, alpha = 0.05)  # power at that budget
checked_rnapower(depth = 20, n = 14, cv = 0.3, effect = 1.5, alpha = 0.05)                  # same n at the deeper depth=20 example
```
At `n = 14, cv = 0.3, effect = 1.5`, the conservative `depth = 2` (20M reads) gives power ~0.29 versus ~0.82 at `depth = 20` — the same replicate count can look adequate or badly underpowered depending on which `depth` was silently assumed, which is why the conversion has to be explicit rather than reusing the vignette's example value.

## Simulation-Based Power -- the Honest Default for Counts

**Goal:** Estimate marginal power and the true realized FDR across the whole expression distribution, accounting for the mean-dispersion trend.

**Approach:** Build (or fit from pilot) a simulation model of counts with a realistic dispersion-mean relationship and DE-effect distribution, simulate many datasets at each candidate sample size, run the intended DE test, and inspect both marginal power **and Actual FDR**. A nominal alpha is only the requested threshold: do not report or select a candidate unless its realized FDR is at or below the pre-specified target (plus an explicitly pre-specified Monte Carlo tolerance). `NaN` Actual FDR is a failed candidate, not a zero FDR.

```r
library(PROPER)
assess_realized_fdr <- function(summary_table, target_fdr = 0.05, tolerance = 0) {
  required <- c("Actual FDR", "Marginal power")
  if (!all(required %in% colnames(summary_table))) {
    stop("summaryPower output is missing Actual FDR or Marginal power.", call. = FALSE)
  }
  actual <- as.numeric(summary_table[, "Actual FDR"])
  out <- data.frame(
    replicates_per_group = if ("SS1" %in% colnames(summary_table)) summary_table[, "SS1"] else seq_along(actual),
    actual_fdr = actual,
    marginal_power = as.numeric(summary_table[, "Marginal power"]),
    accepted = is.finite(actual) & actual <= target_fdr + tolerance
  )
  out$decision <- ifelse(out$accepted, "ACCEPT", "REJECT")
  out
}
sim_opts <- RNAseq.SimOptions.2grp(ngenes = 20000, p.DE = 0.05,
                                    lOD = 'cheung', lBaselineExpr = 'cheung')  # empirical dispersion/expr priors
sims <- runSims(Nreps = c(3, 5, 8, 12), sim.opts = sim_opts, nsims = 50,
                DEmethod = 'edgeR')
powr <- comparePower(sims, alpha.type = 'fdr', alpha.nominal = 0.05,
                     stratify.by = 'expr', delta = log(1.5))          # delta is NATURAL-log lfc in PROPER; marginal power by expression stratum
fdr_gate <- assess_realized_fdr(summaryPower(powr), target_fdr = 0.05, tolerance = 0)
print(fdr_gate)
if (!all(fdr_gate$accepted)) {
  message("REJECTED: do not report or select power for failed candidates. Increase nsims, then refit the pilot dispersion/effect model or revise the DE method before trying again.")
} else {
  plotPower(powr) # grant-ready curve only after every shown candidate passes the realized-FDR gate
}
```

## scRNA-seq Power -- Pseudobulk on Donors, Not Cells

**Goal:** Size a scRNA-seq cross-condition DE study by the quantity that actually sets population power: number of donors, not number of cells.

**Approach:** powsimR is the tool most often named for this, but it is GitHub-only with a compile-required dependency (`bayNorm`) and is not required — population DE power is a donor-level NB power problem, so aggregating (summing) counts per donor into a pseudobulk matrix and running the same edgeR machinery already used elsewhere in this Skill answers it directly, with `n` = number of donors. `examples/scrna_pseudobulk_power.R` runs a full donor x cell simulation end-to-end (edgeR 4.4.2, verified 2026-09-17) and contrasts pseudobulk power against the anti-pattern of testing cells as if they were independent replicates.

The pseudobulk step itself is `DGEList` on the donor-summed matrix (genes x donors, one column per donor) -> `filterByExpr` -> `calcNormFactors` -> `estimateDisp` -> `glmQLFit`/`glmQLFTest`, then BH-adjusted p-values read against the known-DE gene set; the full code is in `examples/scrna_pseudobulk_power.R`. Apply the same realized-FDR gate to each donor-count configuration: a power estimate from a configuration with non-finite or above-target realized FDR is **rejected**, not a candidate design.

Run `Rscript examples/scrna_pseudobulk_power.R` for the full sweep: at 4 donors/group, pseudobulk power stays near 0 regardless of cells per donor (200 vs. 50 makes no difference), while naive cell-level testing looks strong (power > 0.8) at a realized FDR near 0.9 — nine in ten "discoveries" false, because cells are pseudoreplicates, not biological replicates. More donors (not more cells) is the only fix; after that, accept only donor-count configurations whose pseudobulk realized FDR passes the pre-specified gate. A quick closed-form cross-check on the pseudobulk matrix (`checked_rnapower(depth = <post-aggregation depth>, n = <n_donors>, ...)`) applies exactly as in the bulk case once counts are aggregated.

## ATAC-seq / ChIP-seq / Methylation Power -- Same NB Machinery, Per-Region Counts

**Goal:** Size a per-region (peak/CpG) differential accessibility, binding, or methylation study without a separate tool.

**Approach:** `RNASeqPower::rnapower()` and `PROPER` model counts as negative-binomial and are not RNA-seq-specific — RNASeqPower's own vignette is titled "Sample Size for RNA-Seq and **similar** Studies." Feed the same functions per-region counts instead of per-gene counts, using an assay-appropriate `cv` (typically higher than RNA-seq's for ATAC/ChIP due to additional library-prep variability; estimate from pilot peak counts via `DESeq2::estimateDispersions` when possible, as for RNA-seq):

```r
library(RNASeqPower)
# Same call as bulk RNA-seq; depth/cv now describe per-peak or per-CpG coverage and variability.
checked_rnapower(depth = 10, n = 6, cv = 0.5, effect = 1.5, alpha = 0.05)   # per-region power, e.g. ATAC peak
```
For marginal (genome-wide) power across all regions rather than one region, reuse the PROPER simulation block above verbatim with region counts in place of gene counts and an appropriate `lOD`/`lBaselineExpr` fit from pilot data (the built-in `'cheung'` priors are RNA-seq-derived and are a rough stand-in only; a pilot-based mean-dispersion fit is preferred when available).

## Proteomics Power -- Correct for Proteome-Wide Multiplicity

**Goal:** Size a per-protein power calculation (Gaussian after log-transform) so it does not understate the requirement once proteome-wide multiple testing is accounted for.

**Approach:** `pwr::pwr.t.test` per protein is correct as a first step, but plugging in the raw `alpha = 0.05` and stopping there ignores that the study will test every protein on the panel — apply a proteome-wide correction (Bonferroni shown; a simulation-based BH step is the tighter alternative) and report both numbers so they are never confused:

```r
library(pwr)
n_proteins <- 4000
d <- 1.2   # Cohen's d for the target log-abundance fold change

raw <- pwr.t.test(d = d, sig.level = 0.05, power = 0.80, type = 'two.sample')
corrected <- pwr.t.test(d = d, sig.level = 0.05 / n_proteins, power = 0.80, type = 'two.sample')
cat('Raw per-protein alpha=0.05:        n =', ceiling(raw$n), 'per group\n')
cat('Bonferroni-corrected for', n_proteins, 'proteins: n =', ceiling(corrected$n), 'per group\n')
```
Verified on pwr 1.3.0: at `d = 1.2` this prints `n = 12` (raw) vs. `n = 43` (Bonferroni) -- the uncorrected route understates the requirement 3-7x across realistic effect sizes. Missingness (MNAR) is a separate proteomics-specific caveat on top of this (imputation or a missingness-aware test changes the effective `d`), not a substitute for the multiplicity correction.

## Depth vs Replicates -- the Budget Question

For bulk RNA-seq differential expression, sequencing depth shows diminishing returns once it is adequate — Liu, Zhou & White 2014 (*Bioinformatics* 30:301) found the inflection near **~10 million mapped reads** in MCF7 (commonly generalized to a 10-20M band) — whereas adding biological replicates improves power across the whole range. Under a fixed budget, allocate to more biological units before more depth. ATAC/ChIP have their own depth floors (library complexity, peak detection), but the principle holds: biological variance, not read count, limits discovery once depth is adequate.

## CV / Dispersion Guidelines (estimate from pilot when possible)

| Material | Typical biological CV | Source / note |
|----------|----------------------|---------------|
| Cell lines (technical replicates) | 0.1-0.2 | low biological variability |
| Inbred mice | 0.2-0.3 | moderate |
| Primary cells / donor-derived | 0.3-0.4 | donor-dependent |
| Human population samples | 0.3-0.5 | high; Hart 2013 *J Comput Biol* 20:970 default examples |

These are starting points, not substitutes for a pilot estimate; real dispersion is study-specific and a literature CV can be off by a factor of two (estimate via DESeq2/edgeR `estimateDispersions` — see experimental-design/sample-size).

## Per-Method Failure Modes

### Single CV for the whole transcriptome
- **Trigger:** one `cv` plugged into `rnapower()` for all genes.
- **Mechanism:** dispersion varies with mean expression; a single CV mis-states low/high-expressed genes.
- **Symptom:** simulation gives materially different power than the closed form.
- **Fix:** simulation-based power (PROPER/powsimR) from the mean-dispersion trend.

### Observed (post-hoc) power
- **Trigger:** "non-significant, but observed power was 0.3, so add samples."
- **Mechanism:** observed power is a monotone function of the p-value (Hoenig-Heisey 2001).
- **Symptom:** circular reasoning that adds nothing to the CI.
- **Fix:** report effect size + CI; do prospective power for the next study.

### Powering to the expected (or pilot-observed) effect
- **Trigger:** setting the effect to the hoped-for or pilot point estimate.
- **Mechanism:** the pilot estimate is itself noisy; building it in bakes in the winner's curse.
- **Symptom:** chronic underpowering; inflated significant effects (Type-M).
- **Fix:** power to the minimum biologically meaningful effect; propagate pilot variance, not its mean.

### Depth instead of replicates
- **Trigger:** "we will sequence deeper rather than add samples."
- **Mechanism:** past ~10-20M reads, biological variance dominates technical (Liu 2014).
- **Symptom:** deep libraries, still underpowered.
- **Fix:** add biological replicates.

### scRNA-seq power computed on cells
- **Trigger:** "100k cells from 2 patients gives huge power."
- **Mechanism:** population DE power is set by the number of biological samples; cells are pseudoreplicates.
- **Symptom:** power estimate wildly optimistic; results do not replicate.
- **Fix:** pseudobulk over donors (edgeR; see scRNA-seq Power section) — confirmed by simulation: cell-level testing reached 0.87-0.89 realized FDR against a nominal 0.05.

### Proteomics power without multiplicity correction
- **Trigger:** `pwr::pwr.t.test` run per protein at raw `alpha = 0.05` and stopped there.
- **Mechanism:** the study tests every protein on the panel, not one; uncorrected alpha is the wrong operating point.
- **Symptom:** computed n looks small and achievable but the realized proteome-wide FDR is far above nominal.
- **Fix:** Bonferroni (or BH-simulation) correction — see Proteomics Power section; verified 3-7x understatement at realistic panel sizes.

### Computed sample size is not fundable
- **Trigger:** `rnapower()` (or `pwr.t.test`) returns a real but impractically large n (e.g. n=3276 for a 1.05-fold change at 95% power).
- **Mechanism:** the target effect/power combination is genuinely unreachable at any affordable n; the function correctly reports this rather than failing silently.
- **Symptom:** a technically correct number nobody can act on.
- **Fix:** report the power actually achievable at the affordable n, or solve `rnapower()` for the minimum detectable effect at that n (supply `n`/`power`, omit `effect`) instead of chasing an unreachable target.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Power >= 0.80 standard; >= 0.90 for pivotal | convention | tolerable Type-II risk |
| Depth saturates ~10-20M mapped reads for DE | Liu 2014 *Bioinformatics* 30:301 | biological variance then dominates |
| `depth` >= 0.1 x millions of mapped reads, for 85-95% of targets | Hart 2013 (RNASeqPower vignette) | converts a real read budget into `rnapower()`'s `depth` |
| >=6 biological replicates recover most true DE | Schurch 2016 *RNA* 22:839 | n=3 misses many true DE at realistic effects |
| Proteome-wide correction raises n ~3-7x over raw per-protein alpha | verified: pwr 1.3.0, d=1.2, 4000 proteins (n=12 raw vs n=43 Bonferroni) | proteomics needs the same multiplicity framing as per-gene FDR |
| Observed power is a function of the p-value | Hoenig-Heisey 2001 *Am Stat* 55:19 | never use it to interpret a null |
| Type-M exaggeration large in noisy small-n | Gelman-Carlin 2014 *Perspect Psychol Sci* 9:641 | significant effects overstated |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Closed-form and simulation power disagree | single CV vs mean-dispersion trend | use simulation for the reported number |
| "Underpowered (observed power 0.3)" to excuse a null | post-hoc power fallacy | report CI; prospective power only |
| Deep libraries still underpowered | depth over replicates | add biological replicates |
| scRNA-seq power absurdly high | power computed on cells | pseudobulk power over donors (edgeR) |
| Proteomics n looks achievable but panel-wide FDR is inflated | no multiplicity correction on per-protein alpha | Bonferroni/BH-corrected `sig.level` |
| `depth` guessed from a read count with no stated conversion | `depth` is per-gene coverage, not library size | `depth ~= 0.1 x reads(millions)`, majority-of-targets floor |
| Significant effect far larger than literature | winner's curse from underpowering | design analysis (Type-S/Type-M); replicate |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Where did the CV come from?" | estimated from pilot dispersions (DESeq2); literature value used only as a conservative cross-check |
| "Why simulation rather than a formula?" | count power is per-gene; simulation captures the mean-dispersion trend, but its Actual FDR gate must pass before marginal power is usable |
| "Is the study powered?" | only if realized FDR is at or below the pre-specified target (plus stated tolerance) and marginal power >= 0.8 for the minimum meaningful fold change |
| "Why not just sequence deeper?" | depth saturates ~10-20M reads (Liu 2014); replicates added instead |
| "Did you correct for testing 4000 proteins at once?" | yes — Bonferroni/BH-corrected `sig.level`, not the raw per-protein 0.05; both n's reported so they are never confused |
| "Observed power of the null?" | observed power is uninformative (Hoenig-Heisey); CI on the effect reported instead |

## References

- Hart SN, Therneau TM, Zhang Y, Poland GA, Kocher JP. 2013. Calculating sample size estimates for RNA sequencing data. *J Comput Biol* 20:970-978.
- Wu H, Wang C, Wu Z. 2015. PROPER: comprehensive power evaluation for differential expression using RNA-seq. *Bioinformatics* 31:233-241.
- Vieth B, Ziegenhain C, Parekh S, Enard W, Hellmann I. 2017. powsimR: power analysis for bulk and single cell RNA-seq experiments. *Bioinformatics* 33:3486-3488.
- Liu Y, Zhou J, White KP. 2014. RNA-seq differential expression studies: more sequence or more replication? *Bioinformatics* 30:301-304.
- Schurch NJ, Schofield P, Gierliński M, et al. 2016. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22:839-851.
- Hoenig JM, Heisey DM. 2001. The abuse of power: the pervasive fallacy of power calculations for data analysis. *Am Stat* 55:19-24.
- Button KS, Ioannidis JPA, Mokrysz C, Nosek BA, Flint J, Robinson ESJ, Munafò MR. 2013. Power failure: why small sample size undermines the reliability of neuroscience. *Nat Rev Neurosci* 14:365-376.
- Gelman A, Carlin J. 2014. Beyond power calculations: assessing Type S (sign) and Type M (magnitude) errors. *Perspect Psychol Sci* 9:641-651.
- Ioannidis JPA. 2005. Why most published research findings are false. *PLoS Med* 2:e124.

## Related Skills

- sample-size - The inverse problem: minimum replicates for a target power at a target FDR
- randomization-blocking - The experimental unit defines what is replicated; blocking changes error variance
- batch-design - Account for batch/blocking factors in the power model
- differential-expression/deseq2-basics - Estimating dispersions from pilot data for the power model
- single-cell/preprocessing - Pseudobulk model underlying scRNA-seq power
- clinical-biostatistics/power-and-sample-size - Power for regulated clinical-trial endpoints
