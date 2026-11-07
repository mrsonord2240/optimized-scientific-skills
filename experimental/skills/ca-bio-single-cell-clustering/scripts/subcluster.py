#!/usr/bin/env python3
"""Subcluster one cluster: recompute HVGs, scaling, PCA, neighbors and Leiden on the subset only.

Usage:   python subcluster.py in.h5ad out.h5ad cluster_key=leiden cluster=3 [resolution=0.4] [n_pcs=20] [seed=0]
Output:  out.h5ad (the subset) with obs column leiden_sub, backend pinned (flavor='igraph', n_iterations=2, directed=False).
Input:   log-normalized X in the object. Do not pass a scaled or subsetted-to-HVG object.
Needs:   scanpy, igraph. Checked on scanpy 1.12.4 (not yet executed after the recut).
"""
import sys

import scanpy as sc

if len(sys.argv) < 3:
    sys.exit(__doc__)
src, dst = sys.argv[1], sys.argv[2]
opts = dict(cluster_key='leiden', cluster='', resolution='0.4', n_pcs='20', seed='0')
for arg in sys.argv[3:]:
    key, _, value = arg.partition('=')
    if key not in opts:
        sys.exit(f'unknown argument {arg}')
    opts[key] = value
if not opts['cluster']:
    sys.exit('ERROR: pass cluster=<label>')

adata = sc.read_h5ad(src)
if opts['cluster_key'] not in adata.obs:
    sys.exit(f"ERROR: no obs column {opts['cluster_key']}")
mask = adata.obs[opts['cluster_key']].astype(str) == opts['cluster']
if mask.sum() < 50:
    sys.exit(f'ERROR: cluster {opts["cluster"]} has only {int(mask.sum())} cells; too few to subcluster')
sub = adata[mask].copy()
if 'highly_variable' in sub.var:
    sub.var = sub.var.drop(columns=['highly_variable'])
seed = int(opts['seed'])
sc.pp.highly_variable_genes(sub, n_top_genes=2000)
sub = sub[:, sub.var.highly_variable].copy()
sc.pp.scale(sub, max_value=10)
sc.tl.pca(sub, n_comps=30)
sc.pp.neighbors(sub, n_neighbors=15, n_pcs=int(opts['n_pcs']), random_state=seed)
sc.tl.leiden(sub, resolution=float(opts['resolution']), key_added='leiden_sub', flavor='igraph', n_iterations=2,
             directed=False, random_state=seed)
print(sub.obs['leiden_sub'].value_counts().to_string())
sub.write_h5ad(dst)
print(f'wrote {dst}')
