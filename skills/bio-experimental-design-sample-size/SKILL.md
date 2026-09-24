---
name: bio-experimental-design-sample-size
category: Protocol Design
description: Estimates the minimum biological replicates (or cells/donors) for a target power at a target FDR in genomics experiments using ssizeRNA, PROPER, and pilot-data dispersion estimation from DESeq2/edgeR, including a pseudobulk-on-donors route for scRNA-seq cohort sizing. Covers the biological-versus-technical replication distinction (technical replicates do not add degrees of freedom for biological inference), replicate-number-versus-sequencing-depth budgeting, scRNA-seq sample-versus-cell allocation under a pseudobulk model, and the critique that "n=3" is a publication convention rather than a power calculation. Use when budgeting a sequencing experiment, writing the sample-size justification in a grant, estimating replicates from pilot data, allocating a fixed budget between samples and depth, or planning scRNA-seq cohort size. For clinical-trial sample size see clinical-biostatistics/power-and-sample-size; for the power-given-n direction see experimental-design/power-analysis.
tool_type: r
primary_tool: ssizeRNA
license: MIT
author: GPTomics
---

## Version Compatibility

Checked on: ssizeRNA 1.3.3, PROPER 1.38.0, DESeq2 1.46.0, edgeR 4.4.2, pwr 1.3.0 (R 4.4.3 / Bioconductor 3.20).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws an error, introspect the installed package and adapt to the actual API. Notes:
- `ssizeRNA_single()` takes **one mean/dispersion scalar for all genes**; `ssizeRNA_vary()` takes **per-gene vectors**. Passing scalars to `ssizeRNA_vary()` raises `Error in integrate(...) : non-finite function value` on 1.3.3 — use `_single` for a single mean/dispersion and reserve `_vary` for real pilot vectors.
- `res$ssize` from both functions is a **1x3 matrix** `(pi0, ssize, power)`, not a scalar — index it: `res$ssize[, "ssize"]`.
- Both estimators normally return `ssize = NA` silently when no n within `maxN` reaches the target. In ssizeRNA 1.3.3, an extremely low `maxN` can instead raise `argument is of length zero`; treat that exact package edge case as unreachable and emit the same clear `no n <= maxN` message — see "When No n Is Reachable" below.
- `powsimR` is GitHub-only, drifts across versions, and its dependency closure (`bayNorm`) can fail to build against a newer Bioconductor than it was pinned to. If it is not installed or fails to build, use the pseudobulk-on-donors pattern below (`ssizeRNA_vary`/PROPER on donor-level pseudobulk counts) — it needs no extra dependency and is not a lesser substitute, since population DE power is set by donors either way (Squair 2021).
- `PROPER::estParam()` errors with `the condition has length > 1` on R >= 4.0 because a plain `matrix` now has class `c("matrix","array")` and the package's `class(X) %in% c(...)` check was written for R < 4.0. Work around it with `oldClass(X) <- "matrix"` before calling `estParam` (verified on PROPER 1.38.0 / R 4.4.3).

**Setup:**
```r
install.packages(c('BiocManager', 'ssizeRNA', 'pwr'))
BiocManager::install(c('PROPER', 'DESeq2', 'edgeR'))
```

**Worked script:** `examples/sample_size_estimation.R` runs sections 1-4 below end to end on the ssizeRNA/DESeq2/pwr defaults and prints real numbers (not `NA`/`NaN`) — run it first to see the shape of the output before adapting parameters.

# Sample Size for Genomics Experiments

**"How many samples do I need?"** -> Find the smallest number of biological replicates per group that achieves a target marginal power at a target FDR, given the dispersion and effect-size distribution expected for the assay — counting biological units, not measurements.
- R: `ssizeRNA::ssizeRNA_single()` / `ssizeRNA::ssizeRNA_vary()`, `ssizeRNA::check.power()` — FDR-aware NB sample size; pilot dispersions from `DESeq2`/`edgeR`

## The Single Most Important Modern Insight -- The Biological Replicate Is the Unit, and n=3 Is a Convention

Sample size is a count of **biological replicates** — independent experimental units (animals, donors, cultures from independent passages), not measurements. Technical replicates (one library split across lanes, one RNA split into preps) reduce measurement noise but add **no degrees of freedom** for biological inference; averaging them into their biological unit is correct, and selling "n = 3 samples x 3 technical reps = 9" as biological power is a standard error (Blainey, Krzywinski & Altman 2014 *Nat Methods* 11:879). The ubiquitous **"n=3" is a publication convention, not a calculation**: in the 48-vs-48 yeast benchmark, **>=6** biological replicates were needed to recover most true DE genes at realistic effect sizes, and below that the choice of DE tool mattered more than at higher n (Schurch 2016 *RNA* 22:839). Human and primary material, with higher dispersion, need more. For single-cell, the corollary is sharp: population-level DE power is set by the **number of donors**, not the number of cells, because cells are pseudoreplicates — pseudobulk per donor is the correct unit (Squair 2021 *Nat Commun* 12:5692; Murphy & Skene 2022 *Nat Commun* 13:7851).

## Algorithmic Taxonomy

| Approach | Model | Tool | Strength | Fails / costs when |
|----------|-------|------|----------|--------------------|
| Single-parameter NB | one mean/dispersion for all genes | `ssizeRNA::ssizeRNA_single` | quick; transparent; no pilot needed | ignores the mean-dispersion trend |
| FDR-aware NB sample size | NB, per-gene mean/dispersion vectors | `ssizeRNA::ssizeRNA_vary` | controls average power at a true FDR with real gene-to-gene heterogeneity | needs per-gene vectors, not scalars — errors on `integrate()` otherwise |
| Pilot-dispersion simulation | empirical dispersions from pilot | `PROPER` | most defensible; study-specific | requires a pilot dataset |
| Verify a planned n | average power + true FDR at fixed n | `ssizeRNA::check.power` | sanity-checks a budget-driven n | not a search over n; NaN FDR at 0 discoveries needs interpreting |
| scRNA-seq cohort sizing | pseudobulk over donors | `ssizeRNA_vary`/`PROPER` on donor-level pseudobulk counts | counts the right unit (donors); needs no extra dependency beyond DESeq2/edgeR | cell-level sizing is wrong unit |
| Per-feature t-test n | Gaussian (Cohen's d), panel-wide alpha | `pwr::pwr.t.test` with adjusted `sig.level` | proteomics/continuous after transform | wrong for raw counts; wrong at per-feature alpha=0.05 on a multi-feature panel |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Bulk RNA-seq, pilot available | estimate dispersions (DESeq2/edgeR), then `ssizeRNA_vary`/PROPER | study-specific dispersion beats a guess |
| Bulk RNA-seq, no pilot | `ssizeRNA_single` with a literature dispersion, stated as approximate | transparent starting point; `_vary` needs vectors it doesn't have yet |
| Budget already fixed at some n | `check.power` to report achieved power and true FDR | answers "is this n adequate?" |
| scRNA-seq disease vs control | size the number of DONORS on pseudobulk counts (`ssizeRNA_vary`/PROPER) | population power scales with donors, not cells |
| ChIP/ATAC/methylation | NB sample size per region; assay floor as minimum | overdispersed counts; detection floor |
| Proteomics (continuous) | `pwr::pwr.t.test` per protein with panel-wide alpha correction, plus missingness caveat | Gaussian after transform; per-protein alpha=0.05 under-corrects for the panel |
| Have technical replicates | collapse to biological units first | technical reps add no biological df |
| Clinical-trial endpoint | -> clinical-biostatistics/power-and-sample-size | regulated regime |

## FDR-Aware NB Sample Size -- ssizeRNA

**Goal:** Find the minimum biological replicates per group for a target power at a target FDR, accounting for the proportion of DE genes and the mean-dispersion structure.

**No pilot yet -- single mean/dispersion for all genes:** use `ssizeRNA_single`. Use a realistic normalized mean count (order 100-500 for typical bulk RNA-seq depth) — a small toy `mu` can push the answer past a low `maxN` and come back `NA` (see "When No n Is Reachable" below).

```r
library(ssizeRNA)
set.seed(20260918)
res <- ssizeRNA_single(nGenes = 20000, pi0 = 0.95, m = 200,  # m: pseudo sample size for the internal simulation, NOT n per group
                       mu = 200, disp = 0.2,                  # mean count + dispersion (from pilot ideally)
                       fc = 1.5, fdr = 0.05, power = 0.80,
                       maxN = 200)
res$ssize[, "ssize"]                                          # minimum n per group -- res$ssize is a 1x3 matrix (pi0, ssize, power)

# Verify a budget-fixed n: average power and TRUE realized FDR
check.power(nGenes = 20000, pi0 = 0.95, m = 6, mu = 200, disp = 0.2, fc = 1.5, fdr = 0.05, sims = 50)
# A NaN true FDR here means ZERO discoveries at this n, not "FDR unknown" -- see "When No n Is Reachable".
```

**With pilot vectors -- heterogeneous mean/dispersion per gene:** use `ssizeRNA_vary`, fed by the `mu_vec`/`disp_vec` estimated from a pilot below. **Never pass it scalars** — `ssizeRNA_vary(mu = 200, disp = 0.2, ...)` raises `Error in integrate(...) : non-finite function value` on ssizeRNA 1.3.3 regardless of the values chosen; it needs a per-gene vector to integrate over.

```r
library(ssizeRNA)
set.seed(20260918)
res <- ssizeRNA_vary(nGenes = length(mu_vec), pi0 = 0.95,
                     mu = mu_vec, disp = disp_vec,            # VECTORS from the pilot fit below
                     fc = 1.5, fdr = 0.05, power = 0.80, maxN = 200)
res$ssize[, "ssize"]
```

## Pilot Dispersions Drive Honest Sample Size

**Goal:** Replace a guessed CV with a measured dispersion-mean trend from pilot data, as the `mu_vec`/`disp_vec` inputs to `ssizeRNA_vary` above or to PROPER below.

**Approach:** Fit dispersions on the pilot with DESeq2 or edgeR and take per-gene vectors, not a single summary number, into the simulation-based estimator.

```r
library(DESeq2)
set.seed(20260918)
dds <- DESeqDataSetFromMatrix(pilot_counts, pilot_coldata, ~ condition)
dds <- DESeq(dds)
disp_vec <- dispersions(dds)                             # per-gene dispersion estimates
mu_vec   <- rowMeans(counts(dds, normalized = TRUE))     # per-gene normalized mean -- pairs with disp_vec
keep     <- is.finite(disp_vec) & is.finite(mu_vec) & mu_vec > 0
disp_vec <- disp_vec[keep]; mu_vec <- mu_vec[keep]
summary(disp_vec)                                        # use the MEDIAN as a single-number cross-check --
                                                          # the MEAN is pulled up ~2x by a handful of high-dispersion genes
# A literature CV can be off by ~2x; a pilot dispersion is the defensible input.
# Verified on synthetic 2-vs-2 pilot data (planted dispersion 0.35): median DESeq2 estimate 0.30-0.37 (ratio 0.86-1.07x).
```

## Pilot-Data Simulation -- PROPER

**Goal:** Simulate power directly from the pilot's own dispersion/mean distribution rather than a single summary vector fed to `ssizeRNA_vary` — PROPER resamples the actual estimated distribution, not just its per-gene point estimates.

**Approach:** `estParam` characterizes the pilot count matrix; `RNAseq.SimOptions.2grp` builds a simulation config from those estimates plus a target fold change; `runSims` simulates each candidate replicate number; `comparePower` reports power/FDR by replicate number.

Run `scripts/proper_power.R` (args: pilot counts CSV, `reps`, `nsims`, `fc`, `max_genes`, `seed`):

```bash
r.sh scripts/proper_power.R pilot_counts.csv 3,6,10,20 20 1.5   # estParam -> RNAseq.SimOptions.2grp -> runSims -> comparePower
```

The script carries the `oldClass(counts_mat) <- "matrix"` workaround above, sets `delta = log2(fc) - 0.01` (`comparePower` counts a DE gene as a target only if `abs(lfc) > delta`, STRICT, and `runSims` plants every DE gene at exactly `log2(fc)`, so `delta = log2(fc)` leaves ZERO target genes and `power.marginal` is all NaN), and prints `names(powres)` (16 fields; there is no `powerAveraged`, and `$` on a missing name returns NULL silently), `powres$Nreps1`, `rowMeans(powres$power.marginal, na.rm = TRUE)` (`power.marginal` is Nreps x nsims) and `summaryPower(powres)` (PROPER's own table: nominal vs actual FDR, marginal power, avg TD/FD per Nreps).

Read the row where marginal power first reaches 0.80 as the sample size per group; if none does, raise the top `Nreps`. Verified on PROPER 1.38.0 (6v6 synthetic pilot, 2,500 genes, `nsims = 8`): with `delta = log2(1.5) - 0.01` marginal power was 0.006 / 0.060 / 0.22 / 0.61 at n = 3 / 6 / 10 / 20 (actual FDR 0.74 / 0.48 / 0.20 / 0.11), while `delta = log2(1.5)` returned `NaN` for every cell. `power.marginal` can also be `NaN` when `nsims` is very low or a simulation has no true discoveries; raise `nsims` before trusting a cell. Interpret the FDR-aware result as in "When No n Is Reachable".

## scRNA-seq Cohort Sizing -- Pseudobulk on Donors

**Goal:** Size the number of DONORS for population-level differential expression, since cells are pseudoreplicates and cell-level testing inflates false discoveries (see "scRNA-seq sized on cells" below).

**Approach:** Aggregate (sum) each donor's cell-level counts per gene into one pseudobulk sample per donor, then size donors exactly like a bulk RNA-seq study — `ssizeRNA_vary`/PROPER on the pseudobulk dispersions. This needs only DESeq2/edgeR + ssizeRNA/PROPER, already installed for the bulk route, and is the fallback when `powsimR` is not installed (see Version Compatibility).

Run `scripts/pseudobulk_donor_ssize.R` (args: `cell_counts.rds` = named list of one genes x cells matrix per donor, `donor_condition.csv` with columns `donor,condition`, then `fc`, `fdr`, `power`, `maxN`, `seed`):

```bash
r.sh scripts/pseudobulk_donor_ssize.R cell_counts.rds donor_condition.csv 1.5 0.05 0.80 200   # sum per donor -> DESeq2 dispersions -> ssizeRNA_vary
```

It sums (not means) cells per gene per donor, fits DESeq2 on the donor-level pseudobulk, and reports the minimum DONORS per group -- NOT cells.

Verified on synthetic 8-donor pilot (donor dispersion 0.35, 150 cells/donor): pseudobulk median dispersion recovered at 0.32; donor count from `ssizeRNA_vary` = 84 at fc=1.5 (comparable order to the bulk case above -- donors are the same statistical unit as bulk replicates once pseudobulked).

## When No n Is Reachable

`ssizeRNA_single`/`_vary` normally return `ssize = NA` **silently** when no n within `maxN` reaches the target -- this is not a computation failure, it means the search ceiling was too low or the target itself is unreachable. ssizeRNA 1.3.3 has one exception: an extremely low ceiling can raise `argument is of length zero` before constructing that `NA` result. Normalize only that exact package error to the same unreachable outcome, while allowing all other errors through. `check.power` returns `fdr_bh_ave = NaN` when the average number of BH discoveries is zero across simulations -- **a NaN true FDR means zero discoveries, not "FDR unknown."**

1. If `ssize` is `NA`: raise `maxN` (e.g. 30 -> 200 -> 1000) and re-run before concluding anything.
2. If it is still `NA` at a practically fundable `maxN` (a few hundred), report the achieved power at that `maxN` instead of a sample size — do not print `NA` as the answer.
3. As a fallback, sweep the fold change upward at the affordable `n` until power reaches the target, and report that **minimum detectable fold change** instead of a sample size (a 1.2-fold target at 90% power can be unreachable at any fundable n — see Anticipated Reviewer Pushback).

```r
safe_ssize <- function(call) tryCatch(call(), error = function(e) {
  if (identical(conditionMessage(e), "argument is of length zero")) return(NULL)
  stop(e)
})
res <- safe_ssize(function() ssizeRNA_single(..., maxN = maxN))
n <- if (is.null(res)) NA_real_ else res$ssize[, "ssize"]
if (length(n) != 1L || is.na(n)) {
  stop(sprintf("no n <= %d reaches the target; raise maxN or revise fc/dispersion", maxN))
}
```

State the seed and the `sims`/`nsims` count alongside any reported n or power. None of these estimators are deterministic without `set.seed()`: `check.power`'s average power moved between 0.1015 and 0.1098 across unseeded calls at `sims = 20` in prior testing.

## Proteomics Sample Size -- Per-Feature t-test Under Multiplicity

**Goal:** Size a per-protein Gaussian test (Cohen's d, after a variance-stabilizing transform) while controlling the FDR across the whole panel, not just one protein at alpha 0.05.

**Why per-protein alpha=0.05 is wrong:** `pwr.t.test(d = 1.2, sig.level = 0.05, power = 0.80)` returns n=12/group, which looks adequate against the assay-floor table below — but at 5,000 proteins with 10% truly changed and 20% MNAR dropout, n=12 delivers a BH marginal power of only ~0.03, not 0.80. Adjust `sig.level` for the panel size before sizing.

```r
library(pwr)
m_effective <- 5000                                      # proteins actually tested after missingness filtering
pwr.t.test(d = 1.2, sig.level = 0.05 / m_effective, power = 0.80)   # Bonferroni: n ~= 44/group
# Bonferroni is the conservative default here; a simulation-based BH sweep (simulate m_effective
# proteins with the expected proportion changed, run BH, sweep n) gives a tighter n when few
# proteins are truly changed, at the cost of writing the simulation instead of a closed form.
```

The Sample Size by Assay table's Proteomics row (below) is a floor for validating a single already-known marker, not for a proteome-wide discovery panel — for a panel, size from the alpha-adjusted `pwr.t.test` above, not the floor.

## Biological vs Technical Replication

Technical replicates estimate measurement variance; biological replicates estimate the variance that generalizes to the population, and only the latter supports inference about the biology. Average or sum technical replicates into their biological unit before any test. "n = 3 samples x 3 technical reps" is n = 3, not n = 9 (Blainey 2014). This is the sample-size face of the experimental-unit principle (see experimental-design/randomization-blocking).

## Replicates vs Depth Under a Fixed Budget

Once depth is adequate (roughly >=10-20M mapped reads for bulk RNA-seq DE), additional biological replicates buy more power than additional depth (Liu 2014 *Bioinformatics* 30:301). Allocate a fixed budget toward more biological units first. scRNA-seq has an analogous rule at the donor level: more donors beat more cells per donor for population DE, with cells per cell type showing diminishing returns past a few hundred (Squair 2021; Murphy-Skene 2022).

## Sample Size by Assay (floors under favorable conditions, not targets)

| Assay | Practical minimum | For small effects | Source / note |
|-------|-------------------|-------------------|---------------|
| Bulk RNA-seq | 3 (convention) | 6-12 | Schurch 2016 *RNA* 22:839: >=6 recovers most true DE across a realistic FC spectrum -- NOT the same n as 80% power at one fixed FC (see Quantitative Thresholds) |
| scRNA-seq (population DE) | 3 donors | 6+ donors | Squair 2021; donors, not cells, drive power |
| ATAC-seq | 2 | 4-6 | library complexity + peak detection floor |
| ChIP-seq | 2 | 3-4 | IDR reproducibility framework (ENCODE) |
| Proteomics (DIA/TMT) | 3 | 6-10 | single-marker validation floor ONLY; a proteome-wide panel needs the alpha-adjusted `pwr.t.test` above (n~=12-44/group at d=1.2 depending on correction), not this row |
| Methylation (array/WGBS) | 4 | 8-12 | high per-CpG variance |

The "minimum" columns are floors that assume low dispersion and large effects; treat them as the smallest defensible n only after a pilot or literature dispersion supports them. For human-donor cohorts (tumor vs normal, disease vs control), these floors and any failure margin are a starting point, not a substitute for the approved protocol -- see the ethics note under Quantitative Thresholds.

## Per-Method Failure Modes

### Technical replicates counted as biological n
- **Trigger:** "n = 9: 3 samples x 3 technical reps."
- **Mechanism:** technical reps add no biological degrees of freedom (Blainey 2014).
- **Symptom:** over-stated power; results do not generalize.
- **Fix:** collapse technical reps to the biological unit; biological n = 3.

### n=3 by convention
- **Trigger:** choosing 3 because "everyone uses 3."
- **Mechanism:** 3 is a habit, not a calculation; misses many true DE (Schurch 2016).
- **Symptom:** chronic underpowering, irreproducibility.
- **Fix:** size from dispersion + target FDR; expect >=6 for realistic effects, more for human material.

### scRNA-seq sized on cells
- **Trigger:** "100k cells from 2 donors is plenty."
- **Mechanism:** population power scales with donors; cells are pseudoreplicates (Squair 2021).
- **Symptom:** false-discovery-laden DE that does not replicate.
- **Fix:** budget for more donors; size on a pseudobulk model.

### Guessed CV instead of pilot dispersion
- **Trigger:** "human samples are ~0.4, so use 0.4."
- **Mechanism:** real dispersion is study-specific; the guess can be off ~2x.
- **Symptom:** the planned n is wrong by a large factor.
- **Fix:** estimate dispersion from any available pilot (DESeq2/edgeR).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| >=6 biological replicates for bulk RNA-seq DE | Schurch 2016 *RNA* 22:839 | recovers most true DE at realistic effects |
| n=3 is a convention, not a calculation | Schurch 2016 | low power and tool-dependent below 6 |
| Donors, not cells, set scRNA-seq DE power | Squair 2021 *Nat Commun* 12:5692 | cells are pseudoreplicates |
| Technical reps add 0 biological df | Blainey 2014 *Nat Methods* 11:879 | only biological reps generalize |
| Depth saturates ~10-20M reads; add replicates | Liu 2014 *Bioinformatics* 30:301 | biological variance dominates |
| Add 10-20% extra units for failures | common practice | RNA degradation, failed libraries |

**">=6" is not the same question as this Skill's calculation, and the two are ~10x apart in this Skill's own worked example.** Schurch's ">=6" is an empirical *recovery* benchmark averaged over a real spectrum of fold changes (most genes change by less than any single target FC). The `ssizeRNA`/PROPER calculations above answer "what n gives 80% *marginal* power at one *fixed minimum* fold change" -- a stricter question. At the Skill's own example parameters (20,000 genes, disp 0.2, 1.5-fold, FDR 0.05) that calculation returns **n in the mid-40s at mu=200, up to 74 at the toy mu=10 an earlier version of this Skill used**, not 6. Report both numbers when quoting either: Schurch's floor as the absolute minimum for any realistic recovery, and the fixed-FC calculation as the n for a defensible, on-target power claim. Do not average or reconcile them into one number -- they are answers to different questions.

**Ethics/protocol note for human-donor cohorts:** when sizing tumor-vs-normal, disease-vs-control or other human-donor studies, the replicate count -- including the 10-20% failure margin above -- must match what the approved IRB/ethics protocol specifies for that cohort. Recruiting additional donors purely to cover the failure margin has its own consent and recruitment implications; a power calculation justifies the number but does not by itself authorize recruiting it.

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Why this n?" | smallest n reaching marginal power >= 0.8 at FDR 0.05 for the minimum meaningful FC; power curve provided |
| "Where did dispersion come from?" | estimated from pilot (DESeq2); literature value used only as a cross-check |
| "Is n=3 enough?" | no; quote both the Schurch floor and the fixed-FC calculation for this target fold change (see the note under Quantitative Thresholds) — the fixed-FC number is the on-target power claim |
| "Why so many donors for scRNA-seq?" | population DE power scales with donors, not cells (Squair 2021) |
| "Technical replicates?" | collapsed to biological units; they add no biological degrees of freedom |

## References

- Bi R, Liu P. 2016. Sample size calculation while controlling false discovery rate for differential expression analysis with RNA-sequencing experiments. *BMC Bioinformatics* 17:146.
- Schurch NJ, Schofield P, Gierliński M, et al. 2016. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22:839-851.
- Blainey P, Krzywinski M, Altman N. 2014. Points of significance: replication. *Nat Methods* 11:879-880.
- Liu Y, Zhou J, White KP. 2014. RNA-seq differential expression studies: more sequence or more replication? *Bioinformatics* 30:301-304.
- Squair JW, Gautier M, Kathe C, et al. 2021. Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692.
- Murphy AE, Skene NG. 2022. A balanced measure shows superior performance of pseudobulk methods in single-cell RNA-sequencing analysis. *Nat Commun* 13:7851.
- Wu H, Wang C, Wu Z. 2015. PROPER: comprehensive power evaluation for differential expression using RNA-seq. *Bioinformatics* 31:233-241.
- Vieth B, Ziegenhain C, Parekh S, Enard W, Hellmann I. 2017. powsimR: power analysis for bulk and single cell RNA-seq experiments. *Bioinformatics* 33:3486-3488.

## Related Skills

- power-analysis - The power-given-n direction and simulation-based power
- randomization-blocking - The experimental unit defines what is counted as a replicate
- batch-design - Balanced designs assume equal n per group
- differential-expression/deseq2-basics - Estimating pilot dispersions for the sample-size model
- single-cell/preprocessing - Pseudobulk aggregation underlying scRNA-seq cohort sizing
- clinical-biostatistics/power-and-sample-size - Sample size for regulated clinical trials
