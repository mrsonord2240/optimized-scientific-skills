"""Welch tests on log2 feature intensities, with explicit BH adjustment."""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

p = argparse.ArgumentParser(description=__doc__)
for name in ('peaks', 'samples', 'case', 'control', 'output'):
    p.add_argument(name)
a = p.parse_args()
peaks = pd.read_csv(a.peaks).set_index('feature')
meta = pd.read_csv(a.samples, dtype={'sample': str, 'group': str})
if meta['sample'].duplicated().any() or not set(peaks.columns) <= set(meta['sample']):
    raise ValueError('sample IDs must be unique and cover peak columns')
if a.case == a.control:
    raise ValueError('case and control must differ')
values = peaks.apply(pd.to_numeric, errors='raise')
if (values < 0).any().any():
    raise ValueError('negative intensities')
logged = np.log2(values.replace(0, np.nan))
groups = meta.set_index('sample')['group']
case = [c for c in peaks.columns if groups[c] == a.case]
control = [c for c in peaks.columns if groups[c] == a.control]
if len(case) < 3 or len(control) < 3:
    raise ValueError('at least three samples per group required')
rows = []
for feature, row in logged.iterrows():
    x, y = row[case].dropna(), row[control].dropna()
    if len(x) >= 3 and len(y) >= 3:
        rows.append((feature, len(x), len(y), x.mean() - y.mean(), ttest_ind(x, y, equal_var=False).pvalue))
result = pd.DataFrame(rows, columns=['feature', 'n_case', 'n_control', 'log2fc', 'pval'])
result['padj'] = np.nan
valid = result['pval'].notna() & np.isfinite(result['pval'])
if valid.any():
    result.loc[valid, 'padj'] = multipletests(result.loc[valid, 'pval'], method='fdr_bh')[1]
result.to_csv(a.output, index=False)
