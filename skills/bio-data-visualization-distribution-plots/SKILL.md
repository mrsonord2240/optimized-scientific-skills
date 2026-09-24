---
name: bio-data-visualization-distribution-plots
description: Plot per-group distributions of continuous data using boxplots, violins, beeswarms, quasirandom jitter, and ggdist raincloud plots with sample-size honesty (Weissgerber 2015), KDE-bandwidth awareness, and N-aware encoding choices. Use when comparing distributions across a small number of groups — expression per cluster, biomarker per arm, scores per condition — and the bar-of-mean default is misleading.
tool_type: mixed
primary_tool: ggplot2
license: MIT
author: GPTomics
---

## Version Compatibility

The R examples were checked on ggplot2 4.0.3, ggbeeswarm 0.7.3, ggdist 3.3.3, lvplot 0.2.2, and introdataviz 0.0.0.9003. The Python snippets were checked on seaborn 0.13.2, matplotlib 3.11.2, and ptitprince 0.3.1.

- Use `ggdist::stat_halfeye` for R rainclouds. Do not make `gghalves` a required dependency: gghalves 0.1.4 is archived and fails with ggplot2 4.x.
- `ptitprince` offers Scott-style bandwidth selection, not Sheather-Jones. This call was checked on matplotlib 3.11.2; recheck it before moving to matplotlib 3.13 because ptitprince has used deprecated matplotlib orientation paths.
- For seaborn categorical palettes, assign `hue=group` and set `legend=False`; palette without hue is deprecated in seaborn 0.13.

Before adapting a pattern, check `packageVersion('<pkg>')` / `?function_name` in R or `help(module.function)` in Python.

# Distribution Plots

**"Plot the distribution per group"** -> choose an encoding from the per-group nonmissing N, retain raw observations at small N, choose a stated KDE bandwidth when density is warranted, and annotate every plotted group with its N. The default `geom_bar(stat = 'summary')` is misleading: Weissgerber et al. (2015) showed that very different distributions can have identical mean bars.

- R: `ggplot2::geom_boxplot`, `ggplot2::geom_violin`, `ggbeeswarm::geom_quasirandom`, `ggdist::stat_halfeye`, `lvplot::geom_lv`
- Python: `seaborn.boxplot/violinplot/swarmplot/stripplot/boxenplot`, `ptitprince.RainCloud`

## Decision Tree by Per-group N

Count nonmissing observations in the data frame actually supplied to each plot.

| N per group | Recommended | Avoid |
|---|---|---|
| 3-10 | Dot plot or jittered raw points + median bar | Bar of mean |
| 10-29 | Beeswarm or quasirandom + box overlay | KDE/violin shape claims |
| 30-200 | ggdist raincloud or box + jitter | Bare violin with an implicit bandwidth |
| 201-1000 | Letter-value plot or violin with an explicit safe bandwidth | Box alone |
| >1000 | Density (KDE) or histogram + summary statistics | Individual points without rasterisation |

Always annotate N in the caption, x-axis label, or a stratum label. For N < 30, show every point; do not use a violin to imply a reliable density estimate.

## Standard Encodings

### Boxplot + deterministic jitter

```r
ggplot(df, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_point(position = position_jitter(width = 0.2, height = 0, seed = 20260923),
             alpha = 0.5, size = 1) +
  scale_fill_manual(values = c('#0072B2', '#D55E00')) +
  labs(x = NULL, y = 'Expression') + theme_classic()
```

Boxplots show median, IQR, 1.5×IQR whiskers, and outliers, but not density or sample size. Suppress boxplot outliers when raw points are overlaid so a point is not drawn twice. Notches use median +/- 1.58 * IQR / sqrt(N); use them only at N >= 15.

For a text annotation instead of axis labels, return both required aesthetics and count nonmissing values:

```r
stat_summary(geom = 'text', vjust = -0.4,
  fun.data = function(x) data.frame(
    y = max(x, na.rm = TRUE), label = paste0('n=', sum(!is.na(x)))
  ))
```

### Violin with a panel-safe bandwidth

`bw = 'SJ'` is useful when it succeeds, but `stats::bw.SJ()` errors for some tied or constant strata. ggplot2 turns that error into a warning and can blank the entire panel. Choose a single panel-wide selector before plotting: if SJ cannot produce a finite positive bandwidth for every nonmissing group, use `nrd0` for the panel and state the fallback in the caption. A completely constant group has no density to estimate; retain its raw points and do not interpret a violin shape for it.

```r
safe_violin_bw <- function(data, group = 'group', value = 'value') {
  values <- split(data[[value]], data[[group]], drop = TRUE)
  sj_ok <- vapply(values, function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L || length(unique(x)) < 2L) return(FALSE)
    bw <- tryCatch(stats::bw.SJ(x), error = function(e) NA_real_)
    is.finite(bw) && bw > 0
  }, logical(1))
  if (all(sj_ok)) 'SJ' else 'nrd0'
}

bw_used <- safe_violin_bw(df)
ggplot(df, aes(group, value, fill = group)) +
  geom_violin(alpha = 0.7, trim = TRUE, bw = bw_used) +
  geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
  labs(caption = paste('Bandwidth:', bw_used))
```

`nrd0` is Silverman's rule; `nrd` has a larger rule-of-thumb bandwidth (1.06 / 0.9 times `nrd0`) and therefore smooths more, not less. Sheather-Jones remains preferred when the guard permits it. For non-negative or otherwise bounded data, keep `trim = TRUE` or set an explicit axis bound: `trim = FALSE` can draw KDE mass outside the possible range.

### Beeswarm / quasirandom

```r
library(ggbeeswarm)
ggplot(df, aes(group, value, color = group)) +
  geom_quasirandom(method = 'quasirandom', width = 0.3, alpha = 0.7) +
  stat_summary(fun = median, geom = 'crossbar', width = 0.5, color = 'black') +
  scale_color_manual(values = c('#0072B2', '#D55E00'))
```

Quasirandom layouts are deterministic and show every point, making them a strong default from N = 10 through N = 29.

### R raincloud (ggdist, N = 30-200)

Use the standalone, current-stack example rather than an archived gghalves recipe:

```r
source('examples/raincloud_phd.R')
```

It uses `ggdist::stat_halfeye`, a compact boxplot, and seeded raw points, and creates its own deterministic data when run as delivered. For a project frame, replace its `df_med` with columns `group` and `value`, then recompute the local N-label table from that frame.

```python
import ptitprince as pt

# ptitprince has no Sheather-Jones selector; Scott can smooth close modes.
pt.RainCloud(x='group', y='value', data=df, hue='group', palette={
    'Control': '#0072B2', 'Treated': '#D55E00'},
    bw='scott', cut=0, width_viol=0.6, orient='h')
```

Use `cut=0` for bounded data so the displayed KDE does not extend beyond observed values. Treat the Python raincloud as a Scott-bandwidth alternative, not as evidence that the R SJ result will be reproduced.

### Letter-value plot (N >= 201)

```r
library(lvplot)
ggplot(df, aes(group, value, fill = group)) +
  geom_lv(k = 5, alpha = 0.7) +
  scale_fill_manual(values = c('#0072B2', '#D55E00'))
```

Letter-value plots retain tail structure that a standard boxplot collapses. In seaborn, avoid the deprecated palette-only form:

```python
sns.boxenplot(data=df, x='group', y='value', hue='group',
              palette={'Control': '#0072B2', 'Treated': '#D55E00'}, legend=False)
```

### Split violin (exactly two complete, sufficiently large conditions)

Split violins are only valid for an ordered, two-level condition factor. Before plotting, remove clusters that lack either condition or have fewer than 30 nonmissing observations in either condition; do not allow a missing cell to flip the sides in later clusters. The runnable guard and example are in `examples/raincloud_phd.R`.

```r
split_df <- prepare_split_violin_data(df_paired, min_n = 30)
ggplot(split_df, aes(cluster, expression, fill = condition,
                     group = interaction(cluster, condition, lex.order = TRUE))) +
  introdataviz::geom_split_violin(alpha = 0.7, trim = TRUE, bw = safe_violin_bw(
    transform(split_df, group = interaction(cluster, condition)), 'group', 'expression')) +
  geom_boxplot(width = 0.15, position = position_dodge(0.5), outlier.shape = NA) +
  scale_fill_manual(values = c(Control = '#56B4E9', Treatment = '#D55E00'))
```

## Failure Modes and Fixes

| Symptom | Cause | Fix |
|---|---|---|
| Means conceal shape or N | bar of mean +/- SEM | use raw points, quasirandom, raincloud, or box + points according to the table |
| `Computation failed in stat_ydensity()` or an empty panel | `bw.SJ` failed in a tied/constant group | use `safe_violin_bw()` panel-wide; retain raw points for constant strata |
| Violin is too smooth to show a plausible second mode | implicit `nrd0` bandwidth | use the guarded SJ choice and compare to a histogram |
| Density appears below zero for non-negative data | `trim = FALSE` extrapolates KDE | use `trim = TRUE` or enforce an explicit lower bound |
| Box + points shows oversized duplicates | boxplot also draws outliers | `geom_boxplot(outlier.shape = NA)` |
| Notch goes outside hinges | N < 15 | remove notches and show raw points |
| N label errors with missing `y` | text stat returned a label only | use an axis label, or return both `y` and `label` from `fun.data` |
| Split-violin sides swap or a tiny cell vanishes | incomplete grid / fewer than 30 observations | call `prepare_split_violin_data()` and omit that cluster from the split-violin panel |
| R raincloud errors with `argument "layout" is missing` | archived gghalves on ggplot2 4.x | use ggdist `stat_halfeye`; do not pin a current workflow to gghalves |

## Quantitative Thresholds

| Threshold | Value | Source |
|---|---|---|
| Valid notched boxplot | N >= 15 | McGill et al. 1978 guidance |
| Show every raw point | N < 30 | Weissgerber et al. 2015 |
| Raincloud / reliable KDE candidate | N = 30-200 | visualization guidance; inspect histogram and bandwidth |
| Letter-value plot | N >= 201 | Hofmann et al. 2017 |
| Split-violin cell eligibility | both ordered conditions present and N >= 30 each | prevents unstable KDE and side ambiguity |
| Whisker length | 1.5 * IQR | Tukey 1977 |
| Notch length | +/- 1.58 * IQR / sqrt(N) | McGill et al. 1978 |

## References

- Allen M, Poggiali D, Whitaker K, Marshall TR, van Langen J, Kievit RA. 2019. Raincloud plots: a multi-platform tool for robust data visualization. *Wellcome Open Res* 4:63. doi:10.12688/wellcomeopenres.15191.1
- Hofmann H, Wickham H, Kafadar K. 2017. Letter-value plots: boxplots for large data. *J Comput Graph Stat* 26(3):469-477. doi:10.1080/10618600.2017.1305277
- McGill R, Tukey JW, Larsen WA. 1978. Variations of box plots. *Am Stat* 32(1):12-16.
- Sheather SJ, Jones MC. 1991. A reliable data-based bandwidth selection method for kernel density estimation. *J R Stat Soc B* 53(3):683-690.
- Tukey JW. 1977. *Exploratory Data Analysis.* Addison-Wesley.
- Weissgerber TL, Milic NM, Winham SJ, Garovic VD. 2015. Beyond bar and line graphs: time for a new data presentation paradigm. *PLOS Biol* 13(4):e1002128. doi:10.1371/journal.pbio.1002128

## Related Skills

- data-visualization/statistical-annotation - Add p-value brackets to distribution plots
- data-visualization/color-palettes - CVD-safe categorical palettes
- data-visualization/ggplot2-fundamentals - Grammar of graphics base
- single-cell/markers-annotation - Stacked / split violin for scRNA gene-by-cluster
- clinical-biostatistics/effect-measures - Effect size to accompany distribution
