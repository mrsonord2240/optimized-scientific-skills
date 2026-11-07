# Purpose: summarise per-sgRNA efficacy from a JACKS grna results file and list, per gene,
#          the fraction of low-efficacy guides (library re-design candidates).
# Inputs:  <outprefix>_grna_JACKS_results.txt (sgrna, X1, X2) and the guide map used for the run
#          (tab-separated with header, sgRNA and Gene columns). Checked on JACKS 0.2, pandas 2.x.
# Usage:   python efficacy_summary.py grna_results.txt guidemap.txt [--low 0.3] [--sgrna-hdr sgRNA] [--gene-hdr Gene]
#          or: from efficacy_summary import efficacy_summary
import argparse
import pandas as pd


def efficacy_summary(grna_results_path, guidemap_path, low_threshold=0.3,
                     sgrna_hdr='sgRNA', gene_hdr='Gene'):
    '''Summarise per-sgRNA efficacy. The grna file has only sgrna/X1/X2, so genes come from the guide map.'''
    df = pd.read_csv(grna_results_path, sep='\t')
    guidemap = pd.read_csv(guidemap_path, sep='\t', usecols=[sgrna_hdr, gene_hdr])
    df = df.merge(guidemap, left_on='sgrna', right_on=sgrna_hdr, how='left')
    unmapped = df[gene_hdr].isna().sum()
    if unmapped:
        raise ValueError(f'{unmapped} sgRNAs in the results are absent from the guide map; check naming')
    df['low_eff'] = df['X1'] < low_threshold
    summary = {
        'total_guides': len(df),
        'low_efficacy_count': int(df['low_eff'].sum()),
        'low_efficacy_pct': df['low_eff'].mean() * 100,
        'median_efficacy': df['X1'].median(),
        'q25_q75': (df['X1'].quantile(0.25), df['X1'].quantile(0.75)),
    }
    # Per-gene proportion of low-efficacy guides
    by_gene = df.groupby(gene_hdr)['low_eff'].mean().sort_values(ascending=False)
    summary['genes_with_all_low_eff'] = int((by_gene == 1).sum())  # genes where every guide is weak
    return summary, by_gene


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('grna_results')
    ap.add_argument('guidemap')
    ap.add_argument('--low', type=float, default=0.3, help='low-efficacy threshold on X1')
    ap.add_argument('--sgrna-hdr', default='sgRNA')
    ap.add_argument('--gene-hdr', default='Gene')
    ap.add_argument('--out', help='optional TSV path for the per-gene low-efficacy fractions')
    a = ap.parse_args()
    summary, by_gene = efficacy_summary(a.grna_results, a.guidemap, a.low, a.sgrna_hdr, a.gene_hdr)
    for k, v in summary.items():
        print(f'{k}: {v}')
    print(by_gene.head(10).to_string())
    if a.out:
        by_gene.to_csv(a.out, sep='\t', header=['low_eff_fraction'])
