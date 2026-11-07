#!/usr/bin/env python3
"""Preprocess single-cell counts with Scanpy: MAD QC, gene filter, normalize, HVGs, PCA.

Usage:   python preprocess_scanpy.py in out.h5ad [tissue=unknown] [mt_prefix=MT-] [batch_key=] [n_top_genes=2000]
                                     [flavor=seurat_v3] [min_cells=3]
         in = Cell Ranger .h5, a filtered_feature_bc_matrix/ directory, or an .h5ad holding raw integer counts.
         tissue: nuclei|unknown = MAD only; pbmc 8; cardiac 30; hepatic 30; skeletal_muscle 40 (mito hard cap, percent).
         batch_key: .obs column of sample ids; MADs are then per sample and HVGs per batch.
Output:  out.h5ad (layers['counts'], X = log-normalized, HVG flags, X_pca) and <out>.qc.tsv.
Stops (exit nonzero) when: no mito gene matches mt_prefix; a QC MAD is 0; under 80% of barcodes survive; input is not
         integer counts; flavor=seurat_v3 without scikit-misc.
Needs:   scanpy, scipy; scikit-misc for seurat_v3. Checked on scanpy 1.12.4 (not yet executed after the recut).
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import median_abs_deviation

if len(sys.argv) < 3:
    sys.exit(__doc__)
src, dst = sys.argv[1], sys.argv[2]
opts = dict(tissue='unknown', mt_prefix='MT-', batch_key='', n_top_genes='2000', flavor='seurat_v3', min_cells='3')
for arg in sys.argv[3:]:
    key, _, value = arg.partition('=')
    if key not in opts:
        sys.exit(f'unknown argument {arg}')
    opts[key] = value
MITO_CAPS = {'nuclei': None, 'pbmc': 8, 'cardiac': 30, 'hepatic': 30, 'skeletal_muscle': 40, 'unknown': None}
if opts['tissue'] not in MITO_CAPS:
    sys.exit(f"ERROR: tissue must be one of {sorted(MITO_CAPS)}")
mito_cap = MITO_CAPS[opts['tissue']]
batch_key = opts['batch_key'] or None

if os.path.isdir(src):
    adata = sc.read_10x_mtx(src, var_names='gene_ids')  # gene_ids are stable across annotation releases
elif src.endswith('.h5'):
    adata = sc.read_10x_h5(src)
else:
    adata = sc.read_h5ad(src)
adata.var_names_make_unique()
print(f'Raw: {adata.n_obs} cells, {adata.n_vars} genes')

values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X).ravel()
if values.size and not np.allclose(values[:200000], np.round(values[:200000])):
    sys.exit('ERROR: X is not integer counts; preprocessing needs RAW counts.')
if batch_key and batch_key not in adata.obs:
    sys.exit(f'ERROR: batch_key {batch_key} not in .obs')

symbols = adata.var['gene_symbols'] if 'gene_symbols' in adata.var else pd.Series(adata.var_names, index=adata.var_names)
adata.var['mt'] = symbols.astype(str).str.startswith(opts['mt_prefix']).to_numpy()
if not adata.var['mt'].any():
    sys.exit(f"ERROR: no gene symbol starts with {opts['mt_prefix']!r}. Ensembl IDs? Match symbols or set mt_prefix (mouse: mt-).")
adata.var['ribo'] = symbols.astype(str).str.startswith(('RPS', 'RPL')).to_numpy()
adata.var['hb'] = symbols.astype(str).str.contains(r'^HB[ABDEGMQZ]\d*(?!\w)').to_numpy()
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], percent_top=[20], log1p=True, inplace=True)
# inplace defaults to False and returns DataFrames; inplace=True writes .obs/.var

report = []


def is_outlier(values, metric, nmads):
    median = np.median(values)
    mad = median_abs_deviation(values)
    report.append((metric, float(median), float(mad), nmads))
    if mad == 0:
        sys.exit(f'ERROR: MAD collapsed for {metric}; use documented fixed cutoffs instead of MAD filtering')
    return ((values < median - nmads * mad) | (median + nmads * mad < values)).to_numpy()


def flag(metric, nmads):
    column = adata.obs[metric]
    if batch_key is None:
        return is_outlier(column, metric, nmads)
    out = np.zeros(adata.n_obs, dtype=bool)
    for _, positions in adata.obs.groupby(batch_key, observed=True).indices.items():
        out[positions] = is_outlier(column.iloc[positions], metric, nmads)
    return out


outlier = (flag('log1p_total_counts', 5) | flag('log1p_n_genes_by_counts', 5) | flag('pct_counts_in_top_20_genes', 5))
mt_outlier = flag('pct_counts_mt', 3)
if mito_cap is not None:
    mt_outlier = mt_outlier | (adata.obs['pct_counts_mt'].to_numpy() > mito_cap)
keep = ~(outlier | mt_outlier)
survival = float(keep.mean())
pd.DataFrame(report, columns=['metric', 'median', 'mad', 'nmads']).assign(survival=survival).to_csv(
    dst.rsplit('.', 1)[0] + '.qc.tsv', sep='\t', index=False)
if survival < 0.80:
    sys.exit(f'ERROR: only {survival:.1%} of barcodes survive QC; inspect the MADs and use fixed, tissue-aware cutoffs')
adata = adata[keep].copy()
sc.pp.filter_genes(adata, min_cells=int(opts['min_cells']))
print(f'Filtered: {adata.n_obs} cells ({survival:.1%} kept), {adata.n_vars} genes')

adata.layers['counts'] = adata.X.copy()  # raw counts for HVG (seurat_v3), doublets, count models
sc.pp.normalize_total(adata)             # target_sum=None: scale each cell to the dataset median
sc.pp.log1p(adata)

flavor, n_top = opts['flavor'], int(opts['n_top_genes'])
if flavor == 'seurat_v3':
    try:
        import skmisc  # noqa: F401
    except ImportError:
        sys.exit("ERROR: flavor=seurat_v3 needs scikit-misc (pip install scikit-misc)")
    hvg_input = adata
    if batch_key:  # per-batch loess needs genes detected in every batch
        keep_genes = np.ones(adata.n_vars, dtype=bool)
        for _, positions in adata.obs.groupby(batch_key, observed=True).indices.items():
            counts = adata.layers['counts'][positions]
            detected = np.asarray((counts > 0).sum(axis=0)).ravel()
            keep_genes &= detected >= int(opts['min_cells'])
        if not keep_genes.any():
            sys.exit('ERROR: no genes meet min_cells in every batch; lower min_cells or inspect the shallow batch')
        hvg_input = adata[:, keep_genes].copy()
    sc.pp.highly_variable_genes(hvg_input, n_top_genes=n_top, flavor='seurat_v3', layer='counts', batch_key=batch_key)
    adata.var['highly_variable'] = False
    adata.var.loc[hvg_input.var_names, 'highly_variable'] = hvg_input.var['highly_variable'].to_numpy()
else:
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor=flavor, batch_key=batch_key)
print(f"HVGs: {int(adata.var['highly_variable'].sum())}")

sc.pp.pca(adata, mask_var='highly_variable')  # PCA on log-normalized HVGs; scaling is skipped by default
adata.write_h5ad(dst)
print(f'wrote {dst}')
