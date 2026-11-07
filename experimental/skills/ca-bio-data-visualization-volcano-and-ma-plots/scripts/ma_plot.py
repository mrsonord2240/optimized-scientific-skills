"""Reference: matplotlib 3.8+ (run clean on 3.11.2), numpy 1.26+ | Verify API if version differs

ma_plot(res, fdr, ax): res is a pandas DataFrame with baseMean, log2FoldChange, padj (shrunken LFC).
"""
import matplotlib.pyplot as plt
import numpy as np

def ma_plot(res, fdr=0.05, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    sig = (res['padj'] < fdr) & res['padj'].notna()
    ax.scatter(np.log10(res.loc[~sig, 'baseMean']), res.loc[~sig, 'log2FoldChange'],
               c='#999999', s=4, alpha=0.4, rasterized=True)
    ax.scatter(np.log10(res.loc[sig, 'baseMean']), res.loc[sig, 'log2FoldChange'],
               c='#D55E00', s=6, alpha=0.7, rasterized=True)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xlabel(r'$\log_{10}$ mean normalized count')
    ax.set_ylabel(r'$\log_2$ fold change (shrunken)')
    return ax
