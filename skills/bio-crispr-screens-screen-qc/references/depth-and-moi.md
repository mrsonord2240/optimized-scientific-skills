# Sequencing Depth and MOI Verification

## Sequencing Depth Audit

**Goal:** Verify that sequencing depth is sufficient to resolve fold changes at the smallest interesting effect size.

**Approach:** Compute reads/sgRNA per sample and the coefficient of variation (CV) of total reads across samples. Compare against Joung 2017's >100 reads/sgRNA for plasmid QC and >500 for screening, or MAGeCK-VISPR's 300x.

```python
def depth_audit(counts_df):
    '''Verify depth: Joung 2017 recommends >100 reads/sgRNA for plasmid QC and
    >500 for screening; MAGeCK-VISPR uses 300x.'''
    total = counts_df.sum()
    n_sgrnas = len(counts_df)
    depth = total / n_sgrnas
    cv = total.std() / total.mean()
    return pd.DataFrame({'total_reads': total, 'reads_per_sgrna': depth,
                          'depth_grade': np.where(depth < 100, 'FAIL',
                                          np.where(depth < 300, 'CAUTION',
                                          np.where(depth < 500, 'OK', 'EXCELLENT')))}).assign(across_sample_cv=cv)
    # 100 = Joung 2017 plasmid-QC floor; 300 = MAGeCK-VISPR; 500 = Joung 2017 screening
```

**CV interpretation:** CV >0.5 across samples in total reads indicates demultiplexing imbalance or library-pooling error; even if individual samples pass depth thresholds, the relative count is then biased.

## MOI Verification

**Goal:** Confirm that infection occurred at MOI 0.3-0.5 so that ≤1 sgRNA/cell predominates.

**Approach:** From titration plate (control wells with serial-diluted virus), compute infection efficiency, then verify by qPCR of integrated proviral copy number in the screen pool.

| MOI | P(≥1 sgRNA/cell) | P(≥2 sgRNAs/cell) | Cells with 2+ guides as fraction of infected |
|-----|------------------|--------------------|-----------------------------------------------|
| 0.3 | 26% | 4% | 14% |
| 0.5 | 39% | 9% | 23% |
| 1.0 | 63% | 26% | 41% |

**Decision rule:** Always titrate to 0.3. At 0.5, 14-23% of "perturbed" cells carry combinatorial perturbations that confound single-gene scoring. The Poisson math is non-negotiable -- there is no analytical correction for high-MOI confounding.
