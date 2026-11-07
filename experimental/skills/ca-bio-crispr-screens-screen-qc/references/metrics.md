## QC Stage Hierarchy

A pooled screen has six distinct bottlenecks where complexity can collapse. Audit each:

| Stage | Metric | Acceptable threshold | Failure consequence | Detail |
|-------|--------|----------------------|----------------------|--------|
| Plasmid pool | Gini, skew, % zero-count guides | Gini <0.1, skew <2 (Joung 2017 states <10), zero <0.5% | Missing guides cannot be screened; dropout indistinguishable from non-coverage | Below; [failure modes](references/failure-modes.md) |
| Day-0 infection | Library coverage, MOI verification | ≥99% guide detection at 500x cells/sgRNA; MOI 0.3 | Founder effects; polyclonality with high MOI | [depth-and-moi](references/depth-and-moi.md) |
| Selection (puro/blast) | % cells surviving, time-course Gini | 30-40% survival at 5-7 days; Gini drift <0.05 | Selection artifact; fast-growers enriched | Below |
| Endpoint | Replicate correlation, depth | Pearson >=0.8 on log-counts (MAGeCK-VISPR floor), Spearman >0.7-0.8, >500 reads/sgRNA (Joung 2017 screening) | Noise dominates; FDR inflates | Below; [depth-and-moi](references/depth-and-moi.md); [pca-and-composite-score](references/pca-and-composite-score.md) |
| Biological signal | CEGv2 PR-AUC, NEGv1 false-positive rate | PR-AUC >0.7 at FDR 5% (community "passing" convention); CEGv2 enrichment in top 1k | Screen lacks essentiality signal; hits not credible | [essentialome-recovery](references/essentialome-recovery.md) |
| Copy-number artifact | Amplified-region enrichment, sgRNA-cut-count correlation | No correlation between sgRNA off-target count and depletion | False-positive essentiality at amplicons; ERBB2 in HER2+ etc. | [copy-number-bias](references/copy-number-bias.md) |

Each metric below, and in the reference files, quantifies one of these stages.

## Library Representation Metrics

**Goal:** Detect dropout, oversaturation, and library bottlenecks at each sequencing stage.

**Approach:** Compute per-sample zero-count fraction, low-count fraction (<30 reads, the CRISPRcleanR `ccr.NormfoldChanges` default), and percentile-based skew, then track how these change between plasmid -> Day-0 -> endpoint to localize the bottleneck.

```bash
python scripts/library_representation.py screen.count.txt --out library_representation.tsv
```

`scripts/library_representation.py` prints one row per sample (`n_sgrnas_detected`, `pct_zero`, `pct_lowcount`, `median_count`, `p10_count`, `p90_count`, `skew_ratio`) and also exports `library_representation(counts_df)` and `stage_specific_thresholds()` for import. The stage limits it returns (pct_zero_max / skew_max / gini_max: plasmid 0.5 / 2.0 / 0.10, day_0 1.0 / 2.5 / 0.12, endpoint 5.0 / 10.0 / 0.30) follow Joung 2017 (zero-count, skew) and MAGeCK-VISPR (Gini); skew 2.0 is a stricter modern convention, Joung 2017 states <10.

**Interpretation:** Plasmid pool failing Gini <0.1 indicates synthesis or amplification bias; the screen is unfit for use. Endpoint Gini drifting above 0.30 indicates either heavy biological selection (acceptable for strong-phenotype drug screens) or a bottleneck (must be diagnosed). The Day-0 vs plasmid delta isolates whether the issue arose during infection (cloning is unlikely to lose specific guides between extraction and infection -- the change happens in cells).

## Gini Coefficient

**Goal:** Quantify how unevenly reads are distributed across sgRNAs in a single sample.

**Approach:** Sort non-zero counts ascending, compute Gini via the cumulative-fraction formula. Compare against stage-specific thresholds.

```python
def gini(x):
    '''Gini coefficient: 0 = perfect equality, 1 = maximal inequality.
    Uses non-zero counts only; zero-count sgRNAs handled separately by % zero.'''
    x = np.sort(x[x > 0].astype(float))
    if x.size == 0:
        return np.nan
    n = x.size
    cumx = np.cumsum(x)
    return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
```

**Stage-specific thresholds:** only the plasmid Gini <=0.1 is a published cutoff (MAGeCK-VISPR); the remaining grades are operational convention.

| Stage | Excellent | Acceptable | Concerning | Failure |
|-------|-----------|------------|------------|---------|
| Plasmid pool | <0.10 | <0.15 | 0.15-0.20 | >0.20 |
| Day 0 (post-infection) | <0.12 | <0.18 | 0.18-0.25 | >0.25 |
| Endpoint (post-selection) | <0.30 | <0.40 | 0.40-0.55 | >0.55 |

A Gini that climbs from 0.10 (plasmid) to 0.45 (endpoint) is expected when the screen exerts strong selection (drug, lethal-condition). A Gini that climbs to 0.45 without any biological selection (e.g., a control-vs-control timepoint comparison) indicates technical drift.

## Replicate Concordance

**Goal:** Verify that biological/technical replicates agree before testing for between-condition differences.

**Approach:** Compute pairwise Pearson on log10(counts+1) (MAGeCK-VISPR convention) and Spearman ρ on raw rank, between every replicate pair within a condition. Flag any pair below the MAGeCK-VISPR floor of 0.8 Pearson on log-scale.

```python
def replicate_concordance(counts_df, condition_map):
    '''condition_map: {condition_name: [sample_col1, sample_col2, ...]}.'''
    log_counts = np.log10(counts_df + 1)
    rows = []
    for cond, samples in condition_map.items():
        if len(samples) < 2:
            continue
        for i in range(len(samples)):
            for j in range(i+1, len(samples)):
                r_pearson = log_counts[[samples[i], samples[j]]].corr().iloc[0, 1]
                r_spearman = counts_df[[samples[i], samples[j]]].corr(method='spearman').iloc[0, 1]
                rows.append({'condition': cond, 'rep1': samples[i], 'rep2': samples[j],
                             'pearson_log': r_pearson, 'spearman': r_spearman})
    return pd.DataFrame(rows)
```

**Thresholds:** 0.8 Pearson is the MAGeCK-VISPR floor; the stricter grades are operational convention.

| Metric | Excellent | Acceptable | Failure |
|--------|-----------|------------|---------|
| Pearson on log10(counts+1) | >0.95 | >0.85 | <0.80 |
| Spearman on raw ranks | >0.85 | >0.70 | <0.60 |

**When Pearson is high but Spearman is low**, a few outlier sgRNAs are driving correlation (one extreme guide dominates). Inspect the scatterplot; typically caused by PCR jackpotting at a single guide. Hit calling should use a method that ranks (RRA, drugZ) rather than one that fits per-sgRNA fold change directly.


## Input Validation

Before any metric, check the count table and fail with a clear message rather than a traceback:

| Check | Expectation | If it fails |
|-------|-------------|-------------|
| Required columns | an sgRNA identifier column (index) and a `Gene` column, then one numeric column per sample; a default `0..n-1` index means the identifier column was lost | name the missing column; do not guess |
| Dtypes | every sample column numeric | report which column is non-numeric and the first offending value |
| All-zero sample | at least one non-zero count per sample | report the sample as failed at sequencing, and skip (not `nan`-propagate) its Gini and correlation |
| Negative counts | none | reject the file; these are not counts |
| Duplicate sgRNA IDs | none | report the duplicates; MAGeCK-format tables should be unique |

`gini()` above already returns `np.nan` for an all-zero sample rather than raising; keep that guard in
any copy of it.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Plasmid Gini | <0.10 | Li W et al 2015 MAGeCK-VISPR *Genome Biol* 16:281 |
| Plasmid skew ratio (p90/p10) | <10 (Joung 2017); <2 is a stricter modern convention | Joung 2017 *Nat Protoc* 12:828 |
| % zero-count sgRNAs (plasmid) | <0.5% | Joung 2017 *Nat Protoc* 12:828 |
| % zero-count sgRNAs (endpoint) | <1% ideal, <5% tolerated | Li W et al 2015 *Genome Biol* 16:281 |
| Replicate Pearson on log10(counts+1) | >=0.8 (MAGeCK-VISPR floor); >0.95 ideal | Li W et al 2015 MAGeCK-VISPR *Genome Biol* 16:281 |
| Replicate Spearman | >0.70 | Operational convention |
| CEGv2 PR-AUC at FDR 5% | >0.70 passing; >0.85 high quality | Community convention (CEGv2 from Hart 2017 *G3* 7:2719) |
| Reads per sgRNA per sample | >100 plasmid QC and >500 screening (Joung 2017); 300+ (MAGeCK-VISPR) | Joung 2017; Li W et al 2015 |
| Library coverage at infection | 500x cells/sgRNA | Joung 2017; DepMap |
| In-vivo coverage | 50-200x at endpoint | Bottleneck-limited; see [[in-vivo-screens]] |
| MOI at infection | 0.3 strict | Poisson: P(≥2)=4% at 0.3 vs 9% at 0.5 |
| CN-bias Spearman ρ (LFC vs copy number), pre-correction fail | abs(ρ) <0.10 | Operational convention |
| CN-bias Spearman ρ, post-correction target | abs(ρ) <0.05 | Operational convention; the tighter bar after CRISPRcleanR/Chronos |
| CN-bias amplified-vs-diploid LFC gap | > -0.5 | Catches focal amplicons that ρ misses |
