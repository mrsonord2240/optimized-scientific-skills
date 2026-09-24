# PCA, Batch Effects and Composite Score

## PCA and Batch Effect Detection

**Goal:** Visualize whether samples cluster by biology or by batch.

**Approach:** PCA on log10(counts+1); samples should cluster by condition, not by batch/replicate-day/library-lot.

```python
from sklearn.decomposition import PCA

def screen_pca(counts_df, metadata_df, condition_col='condition'):
    '''metadata_df: rows = samples, columns include condition_col, batch (optional).'''
    log_counts = np.log10(counts_df + 1).T  # samples as rows for PCA
    pca = PCA(n_components=3)
    pcs = pca.fit_transform(log_counts)
    out = pd.DataFrame(pcs, columns=['PC1', 'PC2', 'PC3'], index=counts_df.columns)
    out = out.join(metadata_df)
    return out, pca.explained_variance_ratio_
```

**Interpretation:** If PC1 separates batches, see [[batch-correction]]. If PC1 separates conditions cleanly, the screen has interpretable biology. If neither separates anything, the screen has no signal (failed) or is dominated by technical noise.

## Composite DepMap-Style Quality Score

**Goal:** Generate a single quality grade combining all metrics for pipeline gating.

**Approach:** Rescale each metric to a comparable 0-1 direction and average them into a single gate score. Screens scoring <-1 SD are typically excluded from DepMap.

```python
def composite_qc_score(per_sample_qc):
    '''per_sample_qc: one row per sample, joining library_representation() output
    (n_sgrnas_detected, reads_per_sgrna) with gini, pearson_min_replicate, pr_auc
    and n_sgrnas_total.'''
    metrics = {
        'gini_inv': 1 - per_sample_qc['gini'],
        'pearson': per_sample_qc['pearson_min_replicate'],
        'pr_auc': per_sample_qc['pr_auc'],
        'depth_log': np.log10(per_sample_qc['reads_per_sgrna']),
        'detected_frac': per_sample_qc['n_sgrnas_detected'] / per_sample_qc['n_sgrnas_total'],
    }
    return pd.DataFrame(metrics).mean(axis=1)
```

This is a pipeline gate, not a publication metric. DepMap reports `gene effect score quality` (Chronos-derived) separately from screen quality; Pacini 2021 scores the latter with NNMD.
