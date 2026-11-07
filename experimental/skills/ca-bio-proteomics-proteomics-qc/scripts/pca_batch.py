#!/usr/bin/env python
# Purpose: PCA on the normalized survivors and a per-PC test of association with batch.
# Inputs:  normalized_log2 (proteins x samples), sample_info indexed by sample name with a batch column.
# Usage:   python scripts/pca_batch.py proteinGroups.txt sample_annotation.csv [--batch-col batch]
#          or: sys.path.insert(0, "scripts"); from pca_batch import pca_batch_check
# Output:  coords, explained_variance_ratio_, tests (pc, p, status = tested / not_testable).
import argparse
import sys

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import f_oneway

def pca_batch_check(normalized_log2, sample_info, batch_col='batch'):
    # sample_info must be indexed by sample name, e.g. pd.read_csv(...).set_index('sample')
    if set(sample_info.index) != set(normalized_log2.columns):
        raise ValueError('sample_info index must equal the matrix column names (set_index on the sample column)')
    # complete cases only: a row-median fill pulls high-missing (failed) samples to the centre of the PCA
    complete = normalized_log2.dropna(how='any')
    n_samples = complete.shape[1]
    if n_samples < 3 or len(complete) < n_samples:
        raise ValueError(f'too few samples ({n_samples}) or complete proteins ({len(complete)}) for PCA')
    n_pc = min(5, n_samples - 1)
    scaled = StandardScaler().fit_transform(complete.T)
    # 'full' is exact and cheap at QC sizes; the default 'auto' is randomized on a wide matrix (see Common Errors)
    pcs = PCA(n_components=n_pc, svd_solver='full', random_state=0).fit(scaled)
    coords = pd.DataFrame(pcs.transform(scaled), columns=[f'PC{i+1}' for i in range(n_pc)],
                          index=complete.columns).join(sample_info)
    print(f'PCA on {len(complete)} complete proteins of {len(normalized_log2)}')
    tests = []
    for pc in coords.columns[:min(3, n_pc)]:
        groups = [coords[coords[batch_col] == b][pc] for b in coords[batch_col].unique()]
        # a level with one sample makes f_oneway raise 'At least two samples are required; got 1'
        if len(groups) < 2 or min(len(g) for g in groups) < 2:
            print(f'{pc} ~ {batch_col}: NOT TESTABLE, level sizes {[len(g) for g in groups]} '
                  f'(need >=2 levels with >=2 samples each) -- this is not evidence of no batch effect')
            tests.append({'pc': pc, 'p': np.nan, 'status': 'not_testable'})
            continue
        _, p = f_oneway(*groups)
        print(f'{pc} ~ {batch_col}: p={p:.4f}')
        tests.append({'pc': pc, 'p': p, 'status': 'tested'})
    return coords, pcs.explained_variance_ratio_, pd.DataFrame(tests)  # tests.status: tested / not_testable


def main():
    ap = argparse.ArgumentParser(description='PCA and batch check')
    ap.add_argument('protein_groups')
    ap.add_argument('sample_sheet')
    ap.add_argument('--sample-col', default='sample')
    ap.add_argument('--batch-col', default='batch')
    ap.add_argument('--group-col', default='condition', help='used only to refuse a batch fully confounded with it')
    ap.add_argument('--prefix', default='LFQ intensity ')
    ap.add_argument('--out', default='pca_batch_results.tsv', help='sample coordinates plus the per-PC batch test')
    a = ap.parse_args()
    pg = pd.read_csv(a.protein_groups, sep='\t', low_memory=False)
    sheet = pd.read_csv(a.sample_sheet).set_index(a.sample_col)
    absent = [a.prefix + s for s in sheet.index if a.prefix + s not in pg.columns]
    if absent:
        sys.exit(f'REFUSED: columns not in the table: {absent[:3]}')
    if a.group_col in sheet.columns and sheet[a.batch_col].nunique() > 1:
        if (sheet.groupby(a.batch_col)[a.group_col].nunique() == 1).all():
            sys.exit(f'REFUSED: every {a.batch_col} level holds a single {a.group_col}, so batch is confounded with it '
                     'and the design is non-identifiable. Report that and stop; do not correct or test.')
    log2 = np.log2(pd.DataFrame({s: pg[a.prefix + s] for s in sheet.index}).replace(0, np.nan))
    coords, evr, tests = pca_batch_check(log2, sheet, batch_col=a.batch_col)
    print(coords.round(3).to_string())
    print(evr.round(4))
    print(tests.to_string())
    coords = coords.copy()
    for i, v in enumerate(evr):
        coords.loc['__explained_variance__', f'PC{i+1}'] = v
    for _, r in tests.iterrows():
        coords.loc[f'__batch_test_{r.pc}__', 'p'] = r.p
        coords.loc[f'__batch_test_{r.pc}__', 'status'] = r.status
    coords.round(5).to_csv(a.out, sep='\t', index_label='sample')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
