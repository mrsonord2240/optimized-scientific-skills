---
name: bio-experimental-design-multiple-testing
category: Data Analysis
description: Controls error rates across thousands of simultaneous tests in genomics discovery using false-discovery-rate methods (Benjamini-Hochberg 1995; Benjamini-Yekutieli 2001 for arbitrary dependence; Storey q-value with pi0 estimation; local FDR; independent filtering Bourgon 2010; covariate-weighted FDR via IHW Ignatiadis 2016), plus family-wise error control (Bonferroni, Holm) and the GWAS genome-wide threshold. Covers the FDR-versus-FWER choice as the discovery-versus-confirmatory distinction, the dependence assumptions behind BH (PRDS) versus BY, pi0 estimation, and the independent-filtering and false-coverage-rate traps. Use when correcting p-values from genome-wide tests, choosing between BH/BY/q-value/Bonferroni, setting an FDR threshold, applying IHW or independent filtering, or interpreting q-values. For confirmatory trials with few pre-specified endpoints (closed testing, graphical/gatekeeping), see clinical-biostatistics/multiplicity-graphical.
tool_type: mixed
primary_tool: qvalue
goal_approach_exempt: true
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: qvalue 2.34+, IHW 1.30+, R stats (base) p.adjust, statsmodels 0.14+, scipy 1.12+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws an error, introspect the installed package and adapt to the actual API. Note: `statsmodels.stats.multitest.multipletests` defaults to `method='hs'` (Holm-Sidak, an FWER method), NOT Benjamini-Hochberg — always pass `method='fdr_bh'`/`'fdr_by'`/`'bonferroni'`/`'holm'` explicitly.

# Multiple Testing Correction

**"Correct p-values for testing thousands of features"** -> Choose an error rate appropriate to the regime (FDR for discovery, FWER for confirmatory), apply a procedure whose dependence assumptions match the data, and report the adjusted quantity with its interpretation.
- R: `p.adjust(p, method = 'BH')`, `qvalue::qvalue()`, `IHW::ihw()`
- Python: `statsmodels.stats.multitest.multipletests(p, method='fdr_bh')`

## The Single Most Important Modern Insight -- FDR vs FWER Is a Choice About Which Error Matters

The choice between false-discovery-rate and family-wise-error control is not a technicality; it is a statement about which kind of mistake is costly. In **discovery** (20,000 genes, thousands of peaks), tolerating a small, controlled fraction of false positives among the rejections buys enormous power — FDR is the right currency, and Bonferroni would discard nearly every true effect. In **confirmatory** work (a handful of pre-specified endpoints), a single false positive is unacceptable and FWER/closed testing is the standard (that regime lives in clinical-biostatistics/multiplicity-graphical). Two further levers buy back power that plain BH leaves on the table: estimating **pi0** (the proportion of true nulls) turns BH into the more powerful **q-value** (Storey 2002 *J R Stat Soc B* 64:479; Storey & Tibshirani 2003 *PNAS* 100:9440), and weighting hypotheses by an **independent informative covariate** recovers power via **IHW** (Ignatiadis 2016 *Nat Methods* 13:577). The dependence structure matters: BH controls FDR under independence or positive regression dependence (PRDS); under arbitrary or negative dependence use **BY** (Benjamini & Yekutieli 2001 *Ann Stat* 29:1165).

## Algorithmic Taxonomy

| Method | Controls | Dependence assumption | When to use | Tool |
|--------|----------|------------------------|-------------|------|
| Bonferroni | FWER | any | tiny families; confirmatory | `p.adjust(method='bonferroni')` |
| Holm | FWER | any | uniformly beats Bonferroni | `p.adjust(method='holm')` |
| Hochberg / Hommel | FWER | positive dependence | step-up FWER, more power | `p.adjust(method='hochberg'/'hommel')` |
| Benjamini-Hochberg | FDR | independence / PRDS | genome-wide discovery default | `p.adjust(method='BH')` |
| Benjamini-Yekutieli | FDR | arbitrary (incl. negative) | dependence actually expected negative (costs ~50% power) | `p.adjust(method='BY')` |
| Storey q-value | pFDR | independence / weak dependence | many true positives (pi0 << 1); family in the thousands+ | `qvalue::qvalue` |
| Local FDR | posterior null prob | two-groups model | per-feature null probability | `qvalue` ($lfdr); `locfdr` |
| IHW | FDR | covariate independent of null p | informative covariate available; pin `nbins` low | `IHW::ihw` |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Genome-wide DE / peaks, discovery | BH or q-value at FDR 0.05 | controlled false-positive fraction; high power |
| Many true positives expected, family in the thousands+ | q-value (estimates pi0) | more powerful than BH when pi0 << 1; unreliable below ~1000 tests, see failure modes |
| Dependence unknown or positive (e.g. co-regulated modules) | BH (still valid) | mean FDP stayed controlled at 0.039-0.045 across all dependence structures tested; only the run-to-run variance rises |
| Dependence actually expected to be negative | BY | BH's proof doesn't cover negative dependence; BY costs ~50% power here, so confirm the dependence first |
| Informative covariate (mean expr, peak width) | IHW, `nbins` pinned low, run in a retried child process; fall back to BH if it persists | data-driven weights recover power if the covariate is genuinely informative, not just null-independent; see failure modes for the solver crash |
| Per-feature "is this one real?" | local FDR | posterior null probability, not tail average |
| Reporting CIs only on significant hits | FCR-adjusted intervals (below) | naive selected CIs under-cover |
| Small confirmatory gene panel | Bonferroni/Holm | FWER appropriate; power loss acceptable |
| GWAS | genome-wide threshold ~5e-8 | ~1M effective independent tests |
| Confirmatory trial, few endpoints | -> clinical-biostatistics/multiplicity-graphical | closed testing / gatekeeping |
| Applying padj to a finished DE table | -> differential-expression/de-results | method choice here; application there |

## FDR -- Benjamini-Hochberg and the q-value

```r
# Benjamini-Hochberg adjusted p-values (the genome-wide default)
padj <- p.adjust(pvalues, method = 'BH')
sum(padj < 0.05)                                  # discoveries at FDR 5%

# Storey q-value: estimates pi0 (fraction of true nulls) for more power when pi0 << 1
# FAMILY-SIZE FLOOR: the default pi0 smoother needs a family in the thousands to be reliable
# (see "q-value fails on small families" below). Below that, use lambda = 0 or plain BH.
library(qvalue)
qobj <- qvalue(pvalues)
qobj$pi0                                           # estimated proportion of true nulls
q   <- qobj$qvalues                                # min FDR at which each feature is called
lfdr <- qobj$lfdr                                  # local FDR: posterior P(null | statistic)
# lfdr cutoff: Efron's convention is to call features with lfdr < 0.2. The mean lfdr over the called set
# estimates that set's FDR, so lfdr < 0.2 is a much stricter rule than "FDR 20%" (m=18,000: lfdr < 0.2
# called 987 features, mean lfdr 0.043, realized FDP 0.0375, qvalue 2.38.0).
called <- lfdr < 0.2; sum(called); mean(lfdr[called])

# Small-family fallback (tested on all-null families of 20-50 GO terms, qvalue 2.38.0):
qobj_small <- qvalue(pvalues, lambda = 0)          # no spline fit; stable at every size tested,
                                                    # returned the correct pi0 = 1.0 on all-null data
# qvalue(pvalues, pi0.method = 'bootstrap') is also documented as a fallback, but testing here found
# it avoids the *spline* error only to fail on a different one (missing-value error) at a similar
# rate -- it is not a reliable fix at small m. Prefer lambda = 0, or drop to plain BH (no pi0 needed).
```

## Dependence -- BH's Guarantee, BY's Cost

BH's proof covers independence and positive regression dependence (PRDS), not arbitrary or negative
dependence. But across 1,200 replicates tested here (independent, positive-block rho=0.8, and
negative-pair structures), BH's *mean* realized FDP never exceeded nominal -- it held at 0.039-0.045
against a 0.05 target in every structure, including strong negative dependence. What dependence
actually does is widen the *spread*: under positive block dependence (the common co-regulated-gene-
module case) SD(FDP) rose from 0.011 to 0.043 and P(FDP > 0.10) rose from 0.000 to 0.092 -- a real
risk, but invisible in any single realized FDR and not, by itself, evidence that BH has failed.

BY is the right tool when dependence is actually expected to be negative -- BH's proof genuinely does
not cover that case -- but it is expensive: power measured here dropped from 0.605 to 0.288 (53%) for
the same protection. Do not reach for BY merely because dependence is unknown or believed positive;
report the correlation structure instead of paying that cost pre-emptively.

```r
padj_by <- p.adjust(pvalues, method = 'BY')        # valid under any dependence structure; costly
```

## Covariate-Weighted FDR -- IHW

Weight hypotheses by an INDEPENDENT informative covariate (e.g. mean expression), which must be
independent of the p-value under the null. Recovers power vs plain BH. **Do not call `ihw()` directly
in your session**: its LP solver can SEGFAULT (a process crash `tryCatch` cannot catch; 50-75% of runs
at m=18,000, see "IHW segfaults or silently reduces to BH"). Use this wrapper, which runs `ihw()` in a
retried child process and falls back to plain BH:

```r
source('scripts/ihw_safe.R')   # defines ihw_safe(p, covariate, alpha = 0.05, nbins = 5, tries = 3); checked on IHW 1.34.0 / R 4.4.3

# de_table: one row per feature, with columns pvalue and mean_expression
res <- ihw_safe(de_table$pvalue, de_table$mean_expression)
de_table$padj_ihw <- res$padj
res$method; sum(res$padj < 0.05)                   # report which method actually produced the padj
```

`scripts/ihw_safe.R` also runs as a CLI: `Rscript scripts/ihw_safe.R in.csv pvalue mean_expression out.csv [alpha]`.

## Independent Filtering -- Power for Free, If the Filter Is Independent

Filtering out features before testing increases power **only if** the filter statistic is independent of the test statistic under the null (Bourgon, Gentleman & Huber 2010 *PNAS* 107:9546). Overall mean count is independent and is why DESeq2 filters low-count genes automatically; a pre-test on variance or a preliminary t-test is **not** independent and biases the FDR. The DE filtering itself is executed in differential-expression; this skill governs whether a proposed filter is legitimate.

## False Coverage Rate -- CIs on a Selected Set

Reporting a naive (1 - alpha) CI only for the features that passed the FDR filter under-covers: naive
95% CIs on a BH-FDR-0.05-selected set covered the planted truth only 82.4% of the time in testing here
(independently re-derived on a harsher selection: 42.9%). The Benjamini & Yekutieli 2005 construction
(*JASA* 100:71) restores coverage by widening the interval in proportion to the fraction of tests
selected:

```r
q <- 0.05                            # the FDR level used to select
R <- sum(padj < q)                   # number of features selected
m <- length(pvalues)                 # total tests
fcr_level  <- 1 - q * R / m          # adjusted confidence level (Benjamini-Yekutieli 2005)
alpha_fcr  <- 1 - fcr_level
z          <- qnorm(1 - alpha_fcr / 2)
ci_lower   <- estimate[padj < q] - z * se[padj < q]
ci_upper   <- estimate[padj < q] + z * se[padj < q]
```

Tested here: naive coverage 82.4% -> FCR-adjusted 97.6% on the same selected set; the independent
re-derivation went 42.9% -> 95.3%. The width increase is the price of controlling coverage *averaged
over the selected set*, not per feature -- report it as such.

## Python Equivalent (mind the default)

```python
from statsmodels.stats.multitest import multipletests
# DEFAULT method is 'hs' (Holm-Sidak, FWER) -- ALWAYS pass method explicitly.
rej, padj, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')   # Benjamini-Hochberg
rej_by, padj_by, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_by')  # BY
```

## GWAS and the Family-Definition Problem

The genome-wide significance threshold of ~5e-8 is a Bonferroni-style bound for roughly one million effectively independent common-variant tests; Dudbridge & Gusnanto 2008 (*Genet Epidemiol* 32:227) derived ~7.2e-8 for European-ancestry data, near the standard 5e-8. The GWAS test machinery lives in population-genetics/association-testing. More broadly, **what counts as "the family"** of tests is an analyst decision and part of the garden of forking paths: correcting within one contrast, across all contrasts, or across a whole paper are different alpha budgets. Pre-specify the family before seeing results.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| q-value finds many more hits than BH | pi0 << 1 (many true positives) | legitimate if the family is large enough for stable pi0 (thousands+); report pi0 |
| BY far more conservative than BH | BY's baseline power cost (~50%) for arbitrary-dependence protection | justified only if dependence is actually expected negative; otherwise BH remains valid and this is not a reason to switch |
| IHW and BH differ substantially | informative, null-independent covariate -- or IHW silently collapsed to BH (see failure modes) | IHW gain is real if independence and informativeness both hold; check `nbins` and rejections before crediting the gain to IHW |
| Filtering changed the hit count | filter not independent of the test statistic | use a null-independent filter (mean count), not variance/preliminary test |
| Per-feature local FDR high but BH q low | tail-average vs per-feature interpretation | report both; local FDR answers "is THIS one real?" |

## Per-Method Failure Modes

### Bonferroni on a transcriptome
- **Trigger:** Bonferroni across 20,000 genes in a discovery study.
- **Mechanism:** FWER control is far too strict for discovery.
- **Symptom:** almost nothing significant; true effects discarded.
- **Fix:** BH or q-value at a target FDR.

### BH under dependence (mean holds; variance widens)
- **Trigger:** BH on correlated test statistics -- co-regulated gene modules, spatial or pedigree structure -- with dependence unknown or positive.
- **Mechanism:** BH's proof covers independence/PRDS only, but across 1,200 replicates tested (independent, positive-block rho=0.8, negative-pair) the *mean* realized FDP never exceeded nominal (0.039-0.045 vs 0.05). What rises is the *variance*: SD(FDP) 0.011 -> 0.043 and P(FDP>0.10) 0.000 -> 0.092 under positive block dependence.
- **Symptom:** no symptom visible in any one realized FDR -- it is run-to-run instability, not a shift in the mean, and an analyst cannot diagnose it from a single dataset.
- **Fix:** BY guarantees the mean under arbitrary (incl. negative) dependence at a real power cost (measured: 0.605 -> 0.288, 53%). Reserve it for dependence you actually expect to be negative; for unknown or positive dependence, BH remains the better default.

### q-value fails on small families
- **Trigger:** `qvalue(pvalues)` on a family of a few dozen to ~100 tests (a GO term list, a small panel).
- **Mechanism:** the default pi0 estimator fits `smooth.spline` over a lambda grid; on small, especially all-null, families the fit hits missing/infinite values and throws instead of returning.
- **Symptom:** `Error in smooth.spline(lambda, pi0, df = smooth.df) : missing or infinite values in inputs are not allowed`; 38.5% of all-null 20-term replicates and 7.5% of all-null 50-term replicates failed this way in testing (qvalue 2.38.0). Where it did return, pi0-hat ran as low as 0.30-0.58 against a planted truth of 1.0.
- **Fix:** below a family in the thousands, use `qvalue(p, lambda = 0)` (stable at every size tested here) or drop to plain BH, which needs no pi0 estimate and never failed at any size (m = 20 to 20,000). `pi0.method = 'bootstrap'` is documented as an alternative but was not reliable here -- it avoids the spline error only to fail on a different one at a comparable rate.

### IHW segfaults or silently reduces to BH
- **Trigger:** `ihw(pvalue ~ covariate, data = ..., alpha = 0.05)` at default `nbins` ("auto") on a family in the thousands; also observed here at explicit low `nbins` (2-5).
- **Mechanism:** IHW 1.34.0 defaults to `lp_solver = "lpsymphony"`; on this build it crashes the R session outright (SIGSEGV, exit 139) with no R-level error to catch -- `tryCatch` cannot help. Testing found the crash at every `nbins` tried from 2 to the default 12, roughly 50-75% of runs on the same m=18,000 data; pinning `nbins` low reduced the rate but did not eliminate it. Separately, IHW's own automatic bin selection (`nbins <- floor(m/1500)`, capped at 40) collapses to a single bin below m ~ 1500, and at `nbins == 1` the function explicitly reduces to plain BH with uniform weights (its own message: "Only 1 bin; IHW reduces to Benjamini Hochberg") -- the covariate weighting silently does nothing. Above 1 bin, IHW also warns "We recommend that you supply (many) more than 1000 p-values..." whenever any bin holds fewer than 1000 tests.
- **Symptom:** R session dies with no traceback (exit 139); or `ihw()` returns but its rejections/padj exactly match plain BH (the silent nbins==1 collapse).
- **Fix:** use `ihw_safe()` (Covariate-Weighted FDR section; `examples/multiple_testing_correction.R` ships the same pattern). It pins `nbins = 5` -- pinning helps but the crash is stochastic here, not deterministic per input -- and retries in a child process, since retrying in the same R session does nothing (`tryCatch` cannot catch a process crash). If it still won't complete it falls back to plain BH: the power loss is real (this Skill's own comparison: up to 238 vs 204 discoveries, ~17%) but BH always completes and needs no LP solver.

### statsmodels default is not BH
- **Trigger:** `multipletests(p)` expecting Benjamini-Hochberg.
- **Mechanism:** default `method='hs'` (Holm-Sidak, FWER).
- **Symptom:** far fewer significant calls than expected.
- **Fix:** pass `method='fdr_bh'` explicitly.

### Non-independent filtering
- **Trigger:** filter on variance or a preliminary test before the main test.
- **Mechanism:** filter statistic correlated with the test statistic under the null (Bourgon 2010).
- **Symptom:** anti-conservative FDR.
- **Fix:** filter only on a null-independent statistic (overall mean count).

### Selected CIs without FCR adjustment
- **Trigger:** reporting unadjusted CIs only for significant features.
- **Mechanism:** selection induces under-coverage (false coverage rate).
- **Symptom:** intervals too narrow; replication misses.
- **Fix:** Benjamini-Yekutieli (2005) FCR adjustment -- see "False Coverage Rate -- CIs on a Selected Set" above.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| FDR < 0.05 discovery default | Benjamini-Hochberg 1995 *JRSS-B* 57:289 | 5% of calls expected false |
| FDR < 0.10 exploratory | common practice | more leads at higher false fraction |
| q-value uses estimated pi0; needs a family in the thousands | Storey 2002 *JRSS-B* 64:479; this Skill's testing (qvalue 2.38.0) | power gain when pi0 << 1; smooth.spline errors on 36-38% of all-null m=20 replicates |
| BH valid under independence/PRDS; mean FDR holds even under negative dependence tested here | Benjamini-Yekutieli 2001 *Ann Stat* 29:1165; this Skill's testing | switch to BY only for actually-expected negative dependence -- it costs ~50% power |
| GWAS ~5e-8 (7.2e-8 derived) | Dudbridge-Gusnanto 2008 *Genet Epidemiol* 32:227 | ~1M effective tests |
| Filter must be null-independent | Bourgon 2010 *PNAS* 107:9546 | otherwise FDR is biased |
| IHW `nbins="auto"` = floor(m/1500), capped 40; reduces to BH at nbins=1 | IHW 1.34.0 source (`ihw.default`) | m < ~1500 gets 1 bin and no covariate gain; see failure modes for the crash risk above that |

## Common Errors

Symptom-first index into the failure modes above -- causes and fixes are documented once, there.

| Symptom | See |
|---------|-----|
| Almost nothing significant genome-wide | Bonferroni on a transcriptome |
| Run-to-run FDP swings under correlated tests | BH under dependence (mean holds; variance widens) |
| `qvalue()` throws, or pi0 looks implausible, on a small family | q-value fails on small families |
| Far fewer hits than expected in Python | statsmodels default is not BH |
| FDR biased after pre-filtering | Non-independent filtering |
| Replication misses "significant" effects | Selected CIs without FCR adjustment |
| R session dies (exit 139), or IHW rejections exactly match BH | IHW segfaults or silently reduces to BH |

## References

- Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 57:289-300.
- Benjamini Y, Yekutieli D. 2001. The control of the false discovery rate in multiple testing under dependency. *Ann Stat* 29:1165-1188.
- Benjamini Y, Yekutieli D. 2005. False discovery rate-adjusted multiple confidence intervals for selected parameters. *J Am Stat Assoc* 100:71-93.
- Storey JD. 2002. A direct approach to false discovery rates. *J R Stat Soc B* 64:479-498.
- Storey JD, Tibshirani R. 2003. Statistical significance for genomewide studies. *PNAS* 100:9440-9445.
- Efron B. 2008. Microarrays, empirical Bayes and the two-groups model. *Stat Sci* 23:1-22.
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. *PNAS* 107:9546-9551.
- Ignatiadis N, Klaus B, Zaugg JB, Huber W. 2016. Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nat Methods* 13:577-580.
- Dudbridge F, Gusnanto A. 2008. Estimation of significance thresholds for genomewide association scans. *Genet Epidemiol* 32:227-234.

## Related Skills

- power-analysis - The FDR target feeds the power/EDR calculation
- sample-size - Replicate number depends on the FDR threshold chosen here
- batch-design - Surrogate variables change the effective number of tests
- differential-expression/de-results - Where the padj column is applied to a DE table
- population-genetics/association-testing - GWAS genome-wide significance machinery
- pathway-analysis/go-enrichment - Correcting enrichment p-values
- clinical-biostatistics/multiplicity-graphical - Confirmatory FWER / closed testing for trials with few endpoints
