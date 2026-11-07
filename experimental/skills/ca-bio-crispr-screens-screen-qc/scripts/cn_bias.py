#!/usr/bin/env python
"""Copy-number amplicon bias diagnostic (Aguirre 2016 / Munoz 2016): two rules, genome-wide and focal.

Purpose: test whether amplified genes are systematically more depleted than diploid genes.
Inputs:  gene-level LFC table (columns gene, lfc) and per-gene copy number (columns gene, copy_number).
Usage:   python cn_bias.py gene_lfc.tsv copy_number.tsv
Import:  from cn_bias import cn_bias_diagnostic
"""
import argparse
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu


def cn_bias_diagnostic(gene_lfc_df, cn_df):
    '''cn_df: per-gene copy number (from WGS/SNP-array/matched ASCAT).
    Tests whether amplified genes show systematically lower LFC.'''
    merged = gene_lfc_df.merge(cn_df, on='gene')
    bins = pd.qcut(merged['copy_number'], q=5, duplicates='drop')
    bin_lfc = merged.groupby(bins, observed=True)['lfc'].agg(['mean', 'median', 'std', 'count'])
    rho, p = spearmanr(merged['copy_number'], merged['lfc'])
    amplified = merged[merged['copy_number'] > 4]['lfc']
    diploid = merged[merged['copy_number'].between(1.5, 2.5)]['lfc']
    gap, p_gap = np.nan, np.nan
    if len(amplified) >= 3 and len(diploid) >= 3:
        gap = amplified.mean() - diploid.mean()
        p_gap = mannwhitneyu(amplified, diploid, alternative='less').pvalue
    return {'cn_vs_lfc_rho': rho, 'cn_vs_lfc_p': p,
            'n_amplified_genes': len(amplified),
            'amplified_mean_lfc': amplified.mean(),
            'diploid_mean_lfc': diploid.mean(),
            'amplified_vs_diploid_gap': gap,        # negative = amplified genes more depleted
            'p_amplified_more_depleted': p_gap,
            'cn_bias_present': bool((rho < -0.1 and p < 0.01) or (gap < -0.5 and p_gap < 0.01)),
            'per_bin': bin_lfc}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('gene_lfc', help='TSV with columns gene, lfc')
    ap.add_argument('copy_number', help='TSV with columns gene, copy_number')
    args = ap.parse_args()
    res = cn_bias_diagnostic(pd.read_csv(args.gene_lfc, sep='\t'),
                             pd.read_csv(args.copy_number, sep='\t'))
    per_bin = res.pop('per_bin')
    for k, v in res.items():
        print(f'{k}: {v:.4g}' if isinstance(v, float) else f'{k}: {v}')
    print(per_bin.round(3).to_string())
