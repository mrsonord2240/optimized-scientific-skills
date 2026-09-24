"""Non-circular checks for a Scanpy partition after clustering."""

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score


def partition_checks(adata, cluster_key, covariates=('batch', 'sample'),
                     resolution=0.5, n_neighbors=15, n_pcs=30,
                     n_bootstrap=20, seed=0):
    """Print covariate alignment and return mean best-Jaccard stability by cluster.

    `adata` must already contain PCA coordinates and `cluster_key`. Passing these
    checks is necessary but not sufficient for a biological population claim.
    """
    if cluster_key not in adata.obs:
        raise KeyError(f"Missing cluster labels: {cluster_key}")
    if 'X_pca' not in adata.obsm:
        raise KeyError("PCA coordinates are required in adata.obsm['X_pca']")
    if not 0 < resolution:
        raise ValueError('resolution must be positive')

    base = adata.obs[cluster_key].astype(str).to_numpy()
    for covariate in covariates:
        if covariate in adata.obs:
            score = adjusted_rand_score(base, adata.obs[covariate].astype(str))
            print(f'ARI({cluster_key}, {covariate}) = {score:.3f}')

    rng = np.random.default_rng(seed)
    rows = []
    for replicate in range(n_bootstrap):
        keep = rng.choice(adata.n_obs, int(0.8 * adata.n_obs), replace=False)
        boot = adata[keep].copy()
        sc.pp.neighbors(boot, n_neighbors=n_neighbors, n_pcs=n_pcs,
                        random_state=seed + replicate)
        sc.tl.leiden(boot, resolution=resolution, key_added='_bootstrap',
                     flavor='igraph', n_iterations=2, directed=False,
                     random_state=seed + replicate)
        bootstrap = boot.obs['_bootstrap'].astype(str).to_numpy()
        for original in np.unique(base[keep]):
            original_mask = base[keep] == original
            best = max(
                (original_mask & (bootstrap == candidate)).sum()
                / (original_mask | (bootstrap == candidate)).sum()
                for candidate in np.unique(bootstrap))
            rows.append((original, best))

    stability = pd.DataFrame(rows, columns=['cluster', 'jaccard']).groupby('cluster').jaccard.mean()
    print('mean bootstrap Jaccard by cluster:')
    print(stability.sort_values().round(3).to_string())
    return stability
