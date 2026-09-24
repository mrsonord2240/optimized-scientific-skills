#!/usr/bin/env python
"""Second-best sgRNA conservative rule: per-gene LFC of the second-most-extreme guide.

Input:   MAGeCK sgrna_summary.txt (columns Gene, LFC).
Output:  per-gene second_best_lfc and single_guide flag (NaN + True for genes with <2 guides).
Usage:   python second_best_lfc.py mageck.sgrna_summary.txt [--direction neg|pos] [-o second_best.tsv]
         (or: from second_best_lfc import second_best_lfc)
"""
import pandas as pd

def second_best_lfc(sgrna_lfc_df, genes_series, direction='neg'):
    '''Return per-gene LFC of the second-best sgRNA in the direction of interest,
    and flag genes with fewer than 2 sgRNAs. For dropout (direction="neg"),
    second-most-negative LFC. A gene with only one sgRNA has no second guide to
    check at all -- return NaN and single_guide=True for it rather than silently
    falling back to the lone guide's own LFC, which would read as "passing" the
    rule with no corroborating guide involved.'''
    results = []
    for gene in genes_series.unique():
        gene_lfc = sgrna_lfc_df[genes_series == gene].sort_values()
        n = len(gene_lfc)
        if n >= 2:
            second = gene_lfc.iloc[1] if direction == 'neg' else gene_lfc.iloc[-2]
            single = False
        else:
            second = float('nan')
            single = True
        results.append({'gene': gene, 'second_best_lfc': second, 'single_guide': single})
    return pd.DataFrame(results)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('sgrna_summary')
    ap.add_argument('--direction', choices=['neg', 'pos'], default='neg')
    ap.add_argument('-o', '--out', default='second_best_lfc.tsv')
    a = ap.parse_args()
    sg = pd.read_csv(a.sgrna_summary, sep='	')
    res = second_best_lfc(sg['LFC'], sg['Gene'], direction=a.direction)
    res.to_csv(a.out, sep='	', index=False)
    print(f"{len(res)} genes, {int(res['single_guide'].sum())} single-guide; wrote {a.out}")
