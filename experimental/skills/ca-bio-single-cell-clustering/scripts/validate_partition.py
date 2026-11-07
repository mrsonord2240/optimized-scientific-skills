#!/usr/bin/env python3
"""Non-circular checks for a chosen Scanpy partition: covariate alignment and bootstrap stability.

Usage:   python validate_partition.py in.h5ad cluster_key [covariates=batch,sample] [resolution=0.5]
                                      [n_neighbors=15] [n_pcs=30] [n_bootstrap=20] [seed=0]
Output:  prints ARI(cluster, covariate) for each covariate in .obs, then mean best-match Jaccard per cluster;
         writes <in>.stability.tsv. Passing is necessary, not sufficient, for a population claim.
Input:   AnnData with X_pca and the cluster column. resolution must equal the one that made the partition.
Needs:   scanpy, igraph, scikit-learn. Checked on scanpy 1.12.4 (not yet executed after the recut).
"""
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

if len(sys.argv) < 3:
    sys.exit(__doc__)
src, cluster_key = sys.argv[1], sys.argv[2]
opts = dict(covariates='batch,sample', resolution='0.5', n_neighbors='15', n_pcs='30', n_bootstrap='20', seed='0')
for arg in sys.argv[3:]:
    key, _, value = arg.partition('=')
    if key not in opts:
        sys.exit(f'unknown argument {arg}')
    opts[key] = value
covariates = [c for c in opts['covariates'].split(',') if c]
resolution, n_neighbors, n_pcs = float(opts['resolution']), int(opts['n_neighbors']), int(opts['n_pcs'])
n_bootstrap, seed = int(opts['n_bootstrap']), int(opts['seed'])

adata = sc.read_h5ad(src)
if cluster_key not in adata.obs:
    sys.exit(f'ERROR: missing cluster labels: {cluster_key}')
if 'X_pca' not in adata.obsm:
    sys.exit("ERROR: PCA coordinates are required in adata.obsm['X_pca']")
if resolution <= 0:
    sys.exit('ERROR: resolution must be positive')

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
    sc.pp.neighbors(boot, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed + replicate)
    sc.tl.leiden(boot, resolution=resolution, key_added='_bootstrap', flavor='igraph', n_iterations=2,
                 directed=False, random_state=seed + replicate)
    bootstrap = boot.obs['_bootstrap'].astype(str).to_numpy()
    for original in np.unique(base[keep]):
        original_mask = base[keep] == original
        best = max((original_mask & (bootstrap == candidate)).sum() / (original_mask | (bootstrap == candidate)).sum()
                   for candidate in np.unique(bootstrap))
        rows.append((original, best))

stability = pd.DataFrame(rows, columns=['cluster', 'jaccard']).groupby('cluster').jaccard.mean()
print('mean bootstrap Jaccard by cluster:')
print(stability.sort_values().round(3).to_string())
stability.to_csv(src.rsplit('.', 1)[0] + '.stability.tsv', sep='\t')
