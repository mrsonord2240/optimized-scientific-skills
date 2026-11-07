#!/usr/bin/env python3
"""Cluster a preprocessed AnnData with Scanpy: PCA, kNN graph, a Leiden resolution sweep, UMAP.

Usage:   python cluster_scanpy.py in.h5ad out.h5ad [n_pcs=30] [n_neighbors=15] [resolutions=0.2,0.4,0.6,0.8,1.0]
                                  [seed=0] [min_dist=0.3]
Output:  out.h5ad with obs columns leiden_r<res> (backend pinned: flavor='igraph', n_iterations=2, directed=False),
         X_pca and X_umap; out.sweep.tsv with cluster count and size range per resolution.
Input:   normalized, log-transformed X with HVGs chosen. Stops on input that looks like raw counts.
Needs:   scanpy, igraph. Checked on scanpy 1.12.4, igraph 1.0.0 (not yet executed after the recut).
"""
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

if len(sys.argv) < 3:
    sys.exit(__doc__)
src, dst = sys.argv[1], sys.argv[2]
opts = dict(n_pcs='30', n_neighbors='15', resolutions='0.2,0.4,0.6,0.8,1.0', seed='0', min_dist='0.3')
for arg in sys.argv[3:]:
    key, _, value = arg.partition('=')
    if key not in opts:
        sys.exit(f'unknown argument {arg}')
    opts[key] = value
n_pcs, n_neighbors, seed = int(opts['n_pcs']), int(opts['n_neighbors']), int(opts['seed'])
resolutions = [float(r) for r in opts['resolutions'].split(',')]

adata = sc.read_h5ad(src)
x = adata.X
sample = x.data[:100000] if sparse.issparse(x) else np.asarray(x).ravel()[:100000]
if np.all(sample == np.round(sample)) and sample.size and sample.max() > 20:
    sys.exit('ERROR: X looks like raw integer counts; cluster normalized, log-transformed data (preprocessing Skill).')
if 'X_pca' not in adata.obsm:
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
if n_pcs > adata.obsm['X_pca'].shape[1]:
    sys.exit(f"ERROR: n_pcs={n_pcs} exceeds the {adata.obsm['X_pca'].shape[1]} PCs available")

sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed)
rows = []
for res in resolutions:
    key = f'leiden_r{res:g}'
    sc.tl.leiden(adata, resolution=res, key_added=key, flavor='igraph', n_iterations=2, directed=False,
                 random_state=seed)
    sizes = adata.obs[key].value_counts()
    rows.append((res, key, sizes.size, int(sizes.min()), int(sizes.max())))
    print(f'resolution {res:g}: {sizes.size} clusters, sizes {sizes.min()}-{sizes.max()}')
sc.tl.umap(adata, min_dist=float(opts['min_dist']), random_state=seed)

sweep = pd.DataFrame(rows, columns=['resolution', 'column', 'n_clusters', 'min_size', 'max_size'])
sweep.to_csv(dst.rsplit('.', 1)[0] + '.sweep.tsv', sep='\t', index=False)
adata.write_h5ad(dst)
print(f'wrote {dst}')
