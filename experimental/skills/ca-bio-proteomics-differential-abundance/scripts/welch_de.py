#!/usr/bin/env python3
"""Welch t-test + BH differential abundance in Python, for large n only (no variance moderation).

Input:  protein table, raw search-engine TSV (linear `Intensity <sample>` columns, empty or 0 = missing,
        `Reverse` / `Potential contaminant` / `Only identified by site` flags) or a CSV of linear
        intensities (first column = protein ID, one column per sample).
        samples.csv with columns sample, condition. Contrast Test-Reference.
Usage:  python welch_de.py proteins.tsv samples.csv Treatment-Control results.csv [prefix="Intensity "] [min_obs=2]
Output: results CSV (protein, log2fc, pvalue, padj, significant) and <results>_untestable.csv
        (proteins with fewer than min_obs values in a group: undetected, not "not significant").
Checked: not executed. Written against pandas 3.0.5, scipy 1.18.1, statsmodels 0.15.0.
"""
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

FLAGS = ['Reverse', 'Potential contaminant', 'Only identified by site']
LARGE_N = 10  # per-group size below which the lack of variance moderation matters


def read_intensities(path, samples, prefix):
    if re.search(r'\.(tsv|txt|tab)$', path, re.I):
        tab = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
        cols = [prefix + s for s in samples]
        missing = [c for c in cols if c not in tab.columns]
        if missing:
            raise SystemExit(f'Columns not in the table: {missing[:5]}. Set prefix= so prefix + sample name matches a column.')
        keep = pd.Series(True, index=tab.index)
        for flag in FLAGS:
            if flag in tab.columns:
                keep &= tab[flag] != '+'
        id_col = next((c for c in ('Protein IDs', 'Majority protein IDs') if c in tab.columns), tab.columns[0])
        out = tab.loc[keep, cols].apply(pd.to_numeric, errors='coerce')
        out.columns = samples
        out.index = pd.Index(tab.loc[keep, id_col])
        return out
    return pd.read_csv(path, index_col=0)[samples]


def preprocess(intensities):
    log2_data = np.log2(intensities.where(intensities > 0))  # zeros are undetected -> NaN, not -inf
    sample_medians = log2_data.median(axis=0)
    return log2_data - sample_medians + sample_medians.median()  # centre every sample on the global median


def differential_abundance(normalized, case_cols, ctrl_cols, min_obs=2):
    rows, untestable = [], []
    for protein, vals in normalized.iterrows():
        case, ctrl = vals[case_cols].dropna(), vals[ctrl_cols].dropna()
        if len(case) >= min_obs and len(ctrl) >= min_obs:
            _, pval = stats.ttest_ind(case, ctrl, equal_var=False)  # Welch; scipy defaults to Student's
            rows.append({'protein': protein, 'log2fc': case.mean() - ctrl.mean(), 'pvalue': pval})
        else:
            untestable.append({'protein': protein, 'n_case': len(case), 'n_ctrl': len(ctrl)})
    if not rows:
        raise SystemExit(f'No protein has >= {min_obs} non-missing values in both groups; a two-sample test is not possible')
    df = pd.DataFrame(rows)
    df['padj'] = multipletests(df['pvalue'].fillna(1.0), method='fdr_bh')[1]  # default is Holm-Sidak; pass fdr_bh
    df['significant'] = df['padj'] < 0.05
    return df.sort_values('pvalue'), pd.DataFrame(untestable, columns=['protein', 'n_case', 'n_ctrl'])


def main(argv):
    if len(argv) < 4:
        raise SystemExit(__doc__)
    path, samples_path, contrast, out_path = argv[:4]
    opt = dict(a.split('=', 1) for a in argv[4:])
    prefix, min_obs = opt.get('prefix', 'Intensity '), int(opt.get('min_obs', 2))
    samples = pd.read_csv(samples_path)
    case_name, ctrl_name = contrast.split('-', 1)
    case_cols = samples.loc[samples['condition'] == case_name, 'sample'].tolist()
    ctrl_cols = samples.loc[samples['condition'] == ctrl_name, 'sample'].tolist()
    if not case_cols or not ctrl_cols:
        raise SystemExit('Contrast must be Test-Reference with both names in the condition column')
    if min(len(case_cols), len(ctrl_cols)) < LARGE_N:
        print(f'WARNING: fewer than {LARGE_N} samples in a group. No variance moderation: report this as exploratory '
              'or use scripts/limma_de.R.', file=sys.stderr)
    intensities = read_intensities(path, case_cols + ctrl_cols, prefix)
    results, untestable = differential_abundance(preprocess(intensities), case_cols, ctrl_cols, min_obs)
    results.to_csv(out_path, index=False)
    untestable.to_csv(re.sub(r'\.csv$', '_untestable.csv', out_path), index=False)
    print(f'tested {len(results)} | untestable {len(untestable)} | significant (padj < 0.05) {int(results["significant"].sum())}')


if __name__ == '__main__':
    main(sys.argv[1:])
