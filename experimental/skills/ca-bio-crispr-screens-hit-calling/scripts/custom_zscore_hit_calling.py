#!/usr/bin/env python
"""Custom z-score gene-level hit calling (RPM -> per-sgRNA LFC -> gene mean -> z -> BH FDR).

Input:   tab-separated count table (sgRNA id in first column, a Gene column, one column per sample).
Output:  per-gene mean_lfc, std_lfc, n_sgrnas, z, p, fdr sorted by z.
Usage:   python custom_zscore_hit_calling.py counts.tsv --ctrl T0_1,T0_2 --treat T18_1,T18_2              [--ntc-prefix NonTargeting] [-o zscore_hits.tsv]
         (or: from custom_zscore_hit_calling import custom_zscore_hit_calling)
Null comes from non-targeting genes when --ntc-prefix is given, else from all genes (assumes <40% change).
"""
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def custom_zscore_hit_calling(counts_df, ctrl_cols, treat_cols, genes_series, ntc_genes=None):
    '''Z-score gene-level hit calling. If ntc_genes provided, null derived from NTCs only;
    otherwise from all genes (assumes <40% changing).'''
    def rpm(df):
        return df.div(df.sum(axis=0), axis=1) * 1e6
    ctrl_rpm = rpm(counts_df[ctrl_cols])
    treat_rpm = rpm(counts_df[treat_cols])
    lfc_per_sgrna = np.log2((treat_rpm.mean(axis=1) + 1) / (ctrl_rpm.mean(axis=1) + 1))
    gene_lfc = pd.DataFrame({'gene': genes_series, 'lfc': lfc_per_sgrna}).groupby('gene')['lfc'].agg(['mean', 'std', 'count'])
    gene_lfc.columns = ['mean_lfc', 'std_lfc', 'n_sgrnas']
    if ntc_genes is not None:
        null = gene_lfc.loc[gene_lfc.index.isin(ntc_genes), 'mean_lfc']
        null_mean, null_std = null.median(), null.std()
    else:
        null_mean = gene_lfc['mean_lfc'].median()
        null_std = gene_lfc['mean_lfc'].std()
    gene_lfc['z'] = (gene_lfc['mean_lfc'] - null_mean) / null_std
    gene_lfc['p'] = 2 * stats.norm.sf(np.abs(gene_lfc['z']))
    gene_lfc['fdr'] = multipletests(gene_lfc['p'], method='fdr_bh')[1]
    return gene_lfc.sort_values('z')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('counts')
    ap.add_argument('--ctrl', required=True, help='comma-separated control columns')
    ap.add_argument('--treat', required=True, help='comma-separated treatment columns')
    ap.add_argument('--gene-col', default='Gene')
    ap.add_argument('--ntc-prefix', default=None, help='gene names starting with this are the null set')
    ap.add_argument('-o', '--out', default='zscore_hits.tsv')
    a = ap.parse_args()
    df = pd.read_csv(a.counts, sep='	', index_col=0)
    genes = df[a.gene_col]
    ntc = None
    if a.ntc_prefix:
        ntc = set(g for g in genes.unique() if str(g).startswith(a.ntc_prefix))
    res = custom_zscore_hit_calling(df, a.ctrl.split(','), a.treat.split(','), genes, ntc_genes=ntc)
    res.to_csv(a.out, sep='	')
    print(f'{len(res)} genes; {(res["fdr"] < 0.05).sum()} at FDR<0.05; wrote {a.out}')
