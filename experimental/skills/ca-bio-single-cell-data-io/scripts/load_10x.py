#!/usr/bin/env python3
"""Load 10x Cell Ranger output into AnnData and save it as h5ad.

Usage:   python load_10x.py <matrix_dir | file.h5> out.h5ad [var_names=gene_ids] [genome=]
         matrix_dir: raw_feature_bc_matrix/ or filtered_feature_bc_matrix/ (MEX); file.h5: Cell Ranger HDF5.
         var_names: gene_ids (default, Ensembl, stable) or gene_symbols; directory input only. For .h5 the names are
                    symbols and Ensembl IDs are in var['gene_ids'].
         genome: pick one genome from a multi-species .h5.
Output:  out.h5ad (gzip, sparse); prints cells x genes, feature types and whether X is integer.
Keeps non-GEX features (gex_only=False): antibody capture and CRISPR guides stay, split on var['feature_types'].
Needs:   scanpy. Checked on scanpy 1.12.4 (not yet executed after the recut).
"""
import os
import sys

import numpy as np
import scanpy as sc
from scipy import sparse

if len(sys.argv) < 3:
    sys.exit(__doc__)
src, dst = sys.argv[1], sys.argv[2]
opts = dict(var_names='gene_ids', genome='')
for arg in sys.argv[3:]:
    key, _, value = arg.partition('=')
    if key not in opts:
        sys.exit(f'unknown argument {arg}')
    opts[key] = value
if opts['var_names'] not in ('gene_ids', 'gene_symbols'):
    sys.exit('ERROR: var_names must be gene_ids or gene_symbols')

if os.path.isdir(src):
    adata = sc.read_10x_mtx(src, var_names=opts['var_names'], gex_only=False, cache=False)
elif src.endswith('.h5'):
    kwargs = {'genome': opts['genome']} if opts['genome'] else {}
    adata = sc.read_10x_h5(src, gex_only=False, **kwargs)
else:
    sys.exit('ERROR: input must be a matrix directory or an .h5 file')
adata.var_names_make_unique()

print(f'{adata.n_obs} cells x {adata.n_vars} genes')
if 'feature_types' in adata.var:
    print(adata.var['feature_types'].value_counts().to_string())
sample = adata.X.data[:100000] if sparse.issparse(adata.X) else np.asarray(adata.X).ravel()[:100000]
print('X is integer counts:', bool(np.all(sample == np.round(sample))), '| sparse:', sparse.issparse(adata.X))
adata.write_h5ad(dst, compression='gzip')
print(f'wrote {dst}')
