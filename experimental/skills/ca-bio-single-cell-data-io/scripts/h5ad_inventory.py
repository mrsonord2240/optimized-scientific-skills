#!/usr/bin/env python3
"""Print the slot inventory of an h5ad, to diff before and after a conversion.

Usage:   python h5ad_inventory.py data.h5ad
Output:  shape, X sparsity and dtype, layers, obs and var columns with dtypes, obsm, varm, obsp, varp, uns keys,
         and raw (shape and var columns) if present.
Needs:   anndata. Checked on anndata 0.13.3 (not yet executed after the recut).
"""
import sys

import anndata as ad
from scipy import sparse

if len(sys.argv) != 2:
    sys.exit(__doc__)
adata = ad.read_h5ad(sys.argv[1])


def describe(x):
    return f"{'sparse' if sparse.issparse(x) else 'dense'} {x.dtype} {x.shape}"


print(f'shape (cells x genes): {adata.n_obs} x {adata.n_vars}')
print(f'X: {describe(adata.X)}')
print('layers:', {k: describe(v) for k, v in adata.layers.items()} or 'none')
print('obs columns:', {c: str(t) for c, t in adata.obs.dtypes.items()} or 'none')
print('var columns:', {c: str(t) for c, t in adata.var.dtypes.items()} or 'none')
for slot in ('obsm', 'varm', 'obsp', 'varp'):
    print(f'{slot}:', {k: tuple(v.shape) for k, v in getattr(adata, slot).items()} or 'none')
print('uns keys:', sorted(adata.uns.keys()) or 'none')
if adata.raw is not None:
    print(f'raw: {describe(adata.raw.X)}; var columns {list(adata.raw.var.columns)}')
else:
    print('raw: none')
