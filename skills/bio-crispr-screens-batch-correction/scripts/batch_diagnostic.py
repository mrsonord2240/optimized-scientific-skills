#!/usr/bin/env python
"""PCA + variance decomposition: is batch or condition driving each principal component?

Inputs : a tab-separated count table (first column = sgRNA id; non-sample columns such as Gene are
         ignored) and a tab-separated metadata table (first column = sample name, plus batch and
         condition columns).
Outputs: a table of PC1-PC5 (fewer if there are fewer samples) with variance explained and the
         one-way ANOVA F and p for batch and for condition. Batch F > 10x condition F on PC1 means
         batch dominates.
Usage  : python batch_diagnostic.py counts.txt metadata.txt [--batch-col batch] [--condition-col condition]
Import : from batch_diagnostic import batch_diagnostic
Needs  : pandas, numpy, scipy, scikit-learn
"""
import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy import stats

def batch_diagnostic(counts_df, metadata_df, batch_col='batch', condition_col='condition'):
    '''Variance decomposition: report fraction of PC1/PC2 variance attributable to batch vs condition.'''
    log_counts = np.log10(counts_df + 1).T  # samples as rows
    n_pc = min(5, log_counts.shape[0], log_counts.shape[1])
    pca = PCA(n_components=n_pc)
    pcs = pca.fit_transform(log_counts)
    out = pd.DataFrame({
        'PC': range(1, n_pc + 1),
        'var_explained': pca.explained_variance_ratio_,
    })
    pc_df = pd.DataFrame(pcs, columns=[f'PC{i+1}' for i in range(n_pc)], index=counts_df.columns).join(metadata_df)
    for i in range(n_pc):
        pc = pc_df[f'PC{i+1}']
        f_b, p_b = stats.f_oneway(*[pc[pc_df[batch_col] == b] for b in pc_df[batch_col].unique()])
        f_c, p_c = stats.f_oneway(*[pc[pc_df[condition_col] == c] for c in pc_df[condition_col].unique()])
        out.loc[i, 'batch_F'] = f_b
        out.loc[i, 'batch_p'] = p_b
        out.loc[i, 'cond_F'] = f_c
        out.loc[i, 'cond_p'] = p_c
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('counts')
    ap.add_argument('metadata')
    ap.add_argument('--batch-col', default='batch')
    ap.add_argument('--condition-col', default='condition')
    a = ap.parse_args()
    raw = pd.read_csv(a.counts, sep='	', index_col=0)
    meta = pd.read_csv(a.metadata, sep='	', index_col=0)
    samples = [c for c in raw.columns if c in meta.index]
    if not samples:
        raise SystemExit('No count-table column matches a metadata sample name')
    out = batch_diagnostic(raw[samples], meta.loc[samples], a.batch_col, a.condition_col)
    print(out.to_string(index=False))
    r = out.loc[0, 'batch_F'] / out.loc[0, 'cond_F']
    print(f'PC1 batch F / condition F = {r:.2f} ({"batch dominates: correction warranted" if r > 10 else "below the 10x threshold"})')
    sig = [int(pc) for pc, p in zip(out['PC'], out['batch_p']) if p < 0.05]
    print(f'PCs where batch is significant (p < 0.05): {sig if sig else "none"}; batch on a later PC also warrants a look at the PCA plot')


if __name__ == '__main__':
    main()
