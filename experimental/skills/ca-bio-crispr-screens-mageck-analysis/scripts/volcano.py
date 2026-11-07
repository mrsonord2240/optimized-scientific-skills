import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def volcano(gene_summary_path, direction='neg', fdr_threshold=0.05, lfc_threshold=1.0):
    '''direction: "neg" for dropout, "pos" for enrichment.'''
    df = pd.read_csv(gene_summary_path, sep='\t')
    lfc_col = f'{direction}|lfc'
    fdr_col = f'{direction}|fdr'
    fig, ax = plt.subplots(figsize=(9, 7))
    sig = (df[fdr_col] < fdr_threshold) & (np.abs(df[lfc_col]) > lfc_threshold)
    ax.scatter(df.loc[~sig, lfc_col], -np.log10(df.loc[~sig, fdr_col].clip(lower=1e-10)),
                c='lightgray', alpha=0.4, s=10)
    ax.scatter(df.loc[sig, lfc_col], -np.log10(df.loc[sig, fdr_col].clip(lower=1e-10)),
                c='red' if direction == 'neg' else 'blue', alpha=0.7, s=18)
    top = df.loc[sig].nsmallest(15, fdr_col)
    for _, r in top.iterrows():
        ax.annotate(r['id'], (r[lfc_col], -np.log10(max(r[fdr_col], 1e-10))), fontsize=7)
    ax.axhline(-np.log10(fdr_threshold), ls='--', c='black', lw=0.5)
    ax.axvline(0, c='gray', lw=0.5)
    ax.set_xlabel(f'{direction} LFC')
    ax.set_ylabel(f'-log10({direction} FDR)')
    return fig

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--direction', choices=['neg', 'pos'], default='neg')
    args = parser.parse_args()
    volcano(args.input, direction=args.direction).savefig(args.output)
