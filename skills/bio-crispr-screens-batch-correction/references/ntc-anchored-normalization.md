## Control-Sgrna Anchored Normalization

**Goal:** Use non-targeting controls as the per-sample reference so batch shifts cancel.

**Approach:** Scale each sample so its NTC sgRNAs have a constant median. Subsequent fold changes are relative to NTCs in each sample, automatically batch-controlling.

```python
def ntc_anchored_normalize(counts_df, ntc_sgrna_names, target_median=1000):
    '''Scale each sample so its NTC median is target_median. Subsequent LFC is NTC-anchored.'''
    is_ntc = counts_df.index.isin(ntc_sgrna_names)
    ntc_medians = counts_df.loc[is_ntc].median(axis=0)
    scale_factors = target_median / ntc_medians.replace(0, np.nan)
    return counts_df * scale_factors, scale_factors
```

**Critical:** Requires ≥500 NTCs in the library (see [[library-design]]). With fewer, the NTC median is unstable and amplifies noise rather than removing batch; fall back to median normalization.
