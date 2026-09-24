#!/usr/bin/env python
"""Three-method consensus hit calling (MAGeCK RRA + BAGEL2 + drugZ) on the same screen.

Inputs:  MAGeCK gene_summary.txt, BAGEL2 bf output (GENE, BF), drugZ output (GENE, fdr_synth).
Output:  table with per-method hit flags and consensus_count (0-3); prints a
         warning when a pair of hit sets is not enriched for overlap (mismatched comparisons).
Usage:   python consensus_hits.py MAGECK.gene_summary.txt BAGEL.bf.txt DRUGZ.txt [-o consensus.tsv]
         (or: from consensus_hits import consensus_hits)
Thresholds default to SKILL.md's Quantitative Thresholds table (FDR<0.05, BF>6).
For the two-method MAGeCK + BAGEL2 version see ../examples/consensus_hits.py.
"""
import pandas as pd
from scipy.stats import hypergeom

def _check_comparable(merged, hit_cols):
    '''Warn if any pair of method hit-sets shows no statistical enrichment for
    overlap -- the signature of merging results that answer different questions
    (e.g. a real essentiality MAGeCK+BAGEL2 pair merged against a drugZ table from
    an unrelated drug-vs-vehicle screen) rather than genuine method disagreement on
    the same comparison. Verified on real data: matched MAGeCK/BAGEL2 hit sets give
    p=0 (highly enriched overlap); a mismatched drugZ table against either gives
    p=1.0 (no enrichment) -- see the "Consensus across 3 methods is empty" entry in [references/failure-modes.md](references/failure-modes.md).'''
    n = len(merged)
    warnings = []
    for i, col_a in enumerate(hit_cols):
        for col_b in hit_cols[i + 1:]:
            a, b = merged[col_a].fillna(False), merged[col_b].fillna(False)
            k, K, N = int((a & b).sum()), int(a.sum()), int(b.sum())
            if K == 0 or N == 0:
                continue
            p = hypergeom.sf(k - 1, n, K, N)
            if p > 0.05:
                warnings.append(f'{col_a} vs {col_b}: overlap not enriched above chance '
                                 f'(observed={k}, expected~{K * N / n:.1f}, p={p:.3f}) -- '
                                 'check these came from the SAME experimental comparison '
                                 'before trusting consensus.')
    for w in warnings:
        print(f'WARNING: {w}')
    return warnings

def consensus_hits(mageck_path, bagel_path, drugz_path,
                   mageck_fdr_thresh=0.05, bagel_bf_thresh=6, drugz_fdr_thresh=0.05):
    '''Build consensus across MAGeCK / BAGEL2 / drugZ on the same screen.
    Each hit gets a count of supporting methods. Defaults match the Quantitative
    Thresholds table below -- keep this function and examples/consensus_hits.py in
    sync with that table, not with each other.'''
    mageck = pd.read_csv(mageck_path, sep='\t')[['id', 'neg|fdr']].rename(columns={'id': 'gene', 'neg|fdr': 'mageck_neg_fdr'})
    bagel = pd.read_csv(bagel_path, sep='\t')[['GENE', 'BF']].rename(columns={'GENE': 'gene', 'BF': 'bagel_bf'})
    drugz = pd.read_csv(drugz_path, sep='\t')[['GENE', 'fdr_synth']].rename(columns={'GENE': 'gene', 'fdr_synth': 'drugz_synth_fdr'})
    merged = mageck.merge(bagel, on='gene', how='outer').merge(drugz, on='gene', how='outer')
    merged['mageck_hit'] = merged['mageck_neg_fdr'] < mageck_fdr_thresh
    merged['bagel_hit'] = merged['bagel_bf'] > bagel_bf_thresh
    merged['drugz_hit'] = merged['drugz_synth_fdr'] < drugz_fdr_thresh
    _check_comparable(merged, ['mageck_hit', 'bagel_hit', 'drugz_hit'])
    merged['consensus_count'] = (merged[['mageck_hit', 'bagel_hit', 'drugz_hit']].astype(int)).sum(axis=1)
    return merged.sort_values('consensus_count', ascending=False)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mageck'); ap.add_argument('bagel'); ap.add_argument('drugz')
    ap.add_argument('--mageck-fdr', type=float, default=0.05)
    ap.add_argument('--bagel-bf', type=float, default=6)
    ap.add_argument('--drugz-fdr', type=float, default=0.05)
    ap.add_argument('-o', '--out', default='consensus_hits_3method.tsv')
    a = ap.parse_args()
    res = consensus_hits(a.mageck, a.bagel, a.drugz, a.mageck_fdr, a.bagel_bf, a.drugz_fdr)
    res.to_csv(a.out, sep='	', index=False)
    print(res['consensus_count'].value_counts().sort_index(ascending=False).to_string())
    print(f'wrote {a.out}')
