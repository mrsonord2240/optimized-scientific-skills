---
name: bio-data-visualization-statistical-annotation
description: Add p-value brackets, significance asterisks, and effect-size annotations to distribution plots using ggpubr, ggsignif, and statannotations with correct test selection (parametric vs non-parametric vs paired), multiple-testing adjustment, and rendering of negative results. Use when a boxplot/violin/raincloud needs in-figure statistical comparisons between groups.
tool_type: mixed
primary_tool: ggpubr
---

# Statistical Annotation

Render pre-specified group comparisons as brackets with the statistical test that matches the design, adjust the intended comparison family, and report an effect size alongside the p-value. The bracket is only presentation; compute and validate the statistics before drawing it.

## Version Compatibility and Installation

The workflows were executed with R 4.4.3, ggpubr 1.0.0, ggsignif 0.6.4, rstatix 1.1.0, lme4 2.0.1, emmeans 2.0.3, Python 3.12, statannotations 0.7.2, seaborn 0.13.2, scipy 1.18.1, and statsmodels 0.15.0.

```r
install.packages(c('ggplot2', 'ggpubr', 'ggsignif', 'rstatix', 'dplyr', 'tidyr', 'lme4', 'emmeans'))
```

```bash
pip install statannotations seaborn scipy statsmodels pandas matplotlib
```

When versions differ, inspect `packageVersion()` / `?function_name` in R or `pip show` / `help()` in Python. Do not retry an incompatible call unchanged.

## Choose the Test from the Design

Do not trust a plotting helper to infer the design. In ggpubr 1.0.0, `stat_compare_means()` defaults to Wilcoxon for two groups and Kruskal-Wallis for more than two; neither default discovers pairing, nesting, repeated measures, or a planned parametric analysis.

| Design | Recommended test | Runnable route |
|---|---|---|
| 2 unpaired groups, approximately normal | Welch t-test | R `t.test()` / ggpubr `method='t.test'`; Python statannotations `t-test_welch` |
| 2 unpaired groups, skewed or small N | Mann-Whitney / Wilcoxon rank-sum | R `wilcox.test()`; Python `Mann-Whitney` |
| 2 paired groups | Paired t or Wilcoxon signed-rank | `scripts/annotate_paired.R` or `.py`; both validate subject IDs first |
| 3+ groups, approximately normal | ANOVA then Tukey HSD | `scripts/annotate_pairwise.R ... tukey tukey` |
| 3+ groups, non-normal | Kruskal-Wallis then Dunn | `scripts/annotate_pairwise.R ... dunn holm` |
| Nested observations | Linear mixed model, then model-based contrasts | `scripts/annotate_nested.R` |
| Time-course / repeated measures | Repeated-measures model or LMM | `nlme::lme()` / mixed model |
| Two-way factorial | Factorial model with interaction | `aov(y ~ a*b)` or an appropriate mixed model |
| Survival / time-to-event | Log-rank or survival model | `survival::survdiff()`; not a t-test |
| Categorical outcome | Chi-square or Fisher exact | `chisq.test()` / `fisher.test()` |

Normality pre-testing alone is a weak selector: Shapiro-Wilk has little power at small N and flags negligible departures at large N. Use the sampling design, distribution shape, residual diagnostics, robustness needs, and domain assumptions together.

## Adjust the Intended Comparison Family

For all pairwise comparisons among K groups there are K(K-1)/2 tests. At nominal 0.05, the probability of at least one false positive is about 14% for 3 groups, 26% for 4, and 54% for 6 under independence.

- `holm`: default for family-wise error control; uniformly at least as powerful as Bonferroni.
- `bonferroni`: simple, conservative family-wise control.
- `BH` / `fdr`: false-discovery-rate control, commonly used for larger genomic families.

The family is the comparisons actually passed to the correction routine. If only two selected pairs are supplied, the correction uses two tests, not every possible pair. Define the scientific family before filtering results, pass all of it together, and state it in the caption.

## R: Adjust First, Then Annotate

Use the runnable workflow (Wilcoxon is the default mode):

```bash
Rscript scripts/annotate_pairwise.R data.csv annotated.png wilcox holm
```

The script requires `group,value`, computes `pairwise_wilcox_test(..., p.adjust.method='holm')`, independently checks its adjusted values with `p.adjust`, adds bracket positions, and writes both the figure and a `*.results.csv` table. The table includes the test, raw and adjusted p-values, adjustment, family size, `effect_type`, and signed `effect_size`; rank-based modes use rank-biserial r and Tukey mode uses Hedges' g, with positive values meaning the first named group tends higher. Use mode `dunn` after Kruskal-Wallis. Use mode `tukey` with adjustment label `tukey` after ANOVA; it runs `rstatix::tukey_hsd()` and uses Tukey's simultaneous adjustment, not Holm, before `stat_pvalue_manual()`.

Bracket spacing is derived from the family size and the y scale is expanded above the omnibus label. For a specific layout, pass a positive spacing fraction as the fifth argument, for example `... wilcox holm 0.16`. When a full family is still too dense, keep the complete results table but use comparison facets or a compact-letter display for the family-wide result instead of forcing every bracket into one panel; state that the alternate view does not show exact per-pair p-values or effects.

Do not use this pattern for adjusted pairwise labels:

```r
stat_compare_means(comparisons = pairs, p.adjust.method = 'holm')
```

In ggpubr 1.0.0, `p.adjust.method` is not a formal argument of `stat_compare_means()`, and the `comparisons=` route displays unadjusted pairwise p-values. `stat_compare_means()` remains suitable for a single overall Kruskal-Wallis or ANOVA label. Alternatively, `ggpubr::geom_pwc(method='wilcox_test', p.adjust.method='holm', label='p.adj.signif')` performs adjusted pairwise annotation.

For ggsignif, compute the adjusted labels first and pass them manually; `test=` computes raw per-pair tests:

```r
geom_signif(comparisons = pairs,
            annotations = stat_test$p.adj.signif,
            y_position = stat_test$y.position)
```

With `map_signif_level=TRUE`, ggsignif 0.6.4 uses `***`, `**`, `*`, and `NS.` rather than ggpubr's four-star / `ns` convention.

## Python: Supply Adjusted Values Explicitly

```bash
python scripts/annotate_pairwise.py data.csv annotated.png holm
```

The script computes Mann-Whitney p-values and signed rank-biserial r values with SciPy, adjusts the whole family with `statsmodels.stats.multitest.multipletests`, supplies those adjusted values through `Annotator.set_pvalues()` before `annotate()`, and persists `effect_type` plus `effect_size` with the p-value fields.

In statannotations 0.7.2, `comparisons_correction='holm'` and `'BH'` are type-1 corrections: they can append a significance-loss suffix, but the displayed number or stars remain based on raw p-values. Do not present those labels as adjusted. Supplying adjusted p-values explicitly works; the built-in `bonferroni` route also rewrites p-values.

`t-test_ind` is Student's independent t-test. Use `t-test_welch` for Welch's unequal-variance test. R may use an exact Wilcoxon calculation where SciPy selects an asymptotic calculation, so cross-language figures can differ slightly unless the same algorithm is requested explicitly.

## Paired Data: Validate IDs Before Testing

Both rstatix paired tests and statannotations paired tests pair values by row order; they do not accept a subject-ID key. Sorting without checking completeness is insufficient. The paired scripts reject duplicate subject/time rows, require the same ID set in both levels, reshape by `subject_id`, and only then test and plot:

```bash
Rscript scripts/annotate_paired.R paired.csv paired.png holm
python scripts/annotate_paired.py paired.csv paired-python.png holm
```

Inputs require `subject_id,time,value`. The R figure uses `ggpaired(id='subject_id')`; the Python figure draws subject trajectories and annotates the independently computed paired p-value.

## Nested Data: Put the Model-Based Contrast on the Figure

Cells or technical replicates are not independent biological replicates. The nested workflow requires `group,subject_id,value`, verifies each subject belongs to one group, fits `lmer(value ~ group + (1|subject_id))`, derives Holm-adjusted `emmeans` contrasts, and passes those p-values to `stat_pvalue_manual()`:

```bash
Rscript scripts/annotate_nested.R nested.csv nested.png holm
```

`summary(lme4::lmer(...))` alone has no fixed-effect p-value column. Use model-based contrasts such as `emmeans`, or aggregate to one pre-specified summary per biological replicate and test those independent summaries.

## Labels, Exact Values, and Effect Sizes

For ggpubr 1.0.0 and statannotations 0.7.2, the default star bins are inclusive at the boundary: `****` for p <= 1e-4, `***` for p <= 0.001, `**` for p <= 0.01, `*` for p <= 0.05, and `ns` otherwise. ggsignif's default mapping differs as noted above. Treat stars as a display convention, not a substitute for exact adjusted p-values.

ggpubr `label='p.format'` prints formatted values and may print `p < 2e-16` for extremely small p-values. Preserve a results table with the test, effect type and signed effect size, raw p, adjusted p, adjustment method, and comparison-family size. The default R and Python pairwise scripts implement this contract; paired and model-based workflows may require design-specific magnitudes rather than reusing an independent-groups effect.

Report magnitude alongside significance:

- Cohen's d for a normal-location comparison.
- Rank-biserial r, computed from Mann-Whitney U as `2U/(n1*n2) - 1` in the shipped pairwise scripts; this is not Cliff's delta.
- Cliff's delta only when it was actually computed as Cliff's delta.
- A median difference and confidence interval when that is easier to interpret.

## Failure Modes

| Symptom | Cause | Fix |
|---|---|---|
| t-test and rank test disagree on skewed small-N data | `method='t.test'` was chosen without checking assumptions | Choose deliberately from design, shape, diagnostics, and robustness needs |
| Raw pairwise stars survive but adjusted p-values do not | Annotation helper displayed raw p | Compute the family once, adjust, then supply `p.adj` / adjusted values |
| Paired plot lines look right but p-value changes after row shuffle | Test paired by row order | Validate identical ID sets and reshape by ID before testing |
| Cell-level p is extreme but patient-level/model p is null | Pseudoreplication | LMM/cluster-aware method or one pre-specified summary per replicate |
| Bracket pyramid omits null comparisons | Selective display | Show all pre-specified comparisons or name the tested subset in the caption |
| Tiny effect has an extreme p at large N | Significance is not magnitude | Report effect size and uncertainty |
| Asterisks have no recoverable value | Exact results were discarded | Save a results table or show formatted adjusted p-values |

## References

- Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate. *J R Stat Soc B* 57:289-300.
- Dunn OJ. 1964. Multiple comparisons using rank sums. *Technometrics* 6:241-252.
- Holm S. 1979. A simple sequentially rejective multiple test procedure. *Scand J Stat* 6:65-70.
- Wasserstein RL, Lazar NA. 2016. The ASA's statement on p-values. *Am Stat* 70:129-133.

## Related Skills

- data-visualization/distribution-plots - underlying box/violin/raincloud plots
- clinical-biostatistics/categorical-tests - categorical outcomes
- clinical-biostatistics/effect-measures - effect sizes and uncertainty
- experimental-design/multiple-testing - choosing and documenting error control
