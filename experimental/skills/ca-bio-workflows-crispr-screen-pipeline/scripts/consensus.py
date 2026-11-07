#!/usr/bin/env python3
"""Tier consensus across hit-calling methods: one row per gene, one hit flag per method, a tier.

Usage: python consensus.py <out.csv> [mageck=<x.gene_summary.txt>] [bagel=<bf.txt>] [drugz=<drugz_out.txt>]
                           [fdr=0.05] [bf=6]
  Give at least two of mageck=, bagel=, drugz=. Every input must be a depletion (negative selection) result.
  hit rules: MAGeCK neg|fdr < fdr; BAGEL2 BF > bf; drugZ fdr_synth < fdr
Output: out.csv with columns gene, <method>_hit, n_hits, tier (Tier-1 = all three, Tier-2 = two, Tier-3 = one).
  With only two methods the best tier is Tier-2. Counts per tier are printed.
Checked: pandas 3.0.5.
"""
import sys

import pandas as pd


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    out = argv[0]
    opt = {'mageck': '', 'bagel': '', 'drugz': '', 'fdr': '0.05', 'bf': '6'}
    for item in argv[1:]:
        key, _, value = item.partition('=')
        if key not in opt:
            sys.exit(f'unknown option: {item}')
        opt[key] = value
    fdr, bf = float(opt['fdr']), float(opt['bf'])
    frames = []
    if opt['mageck']:
        t = pd.read_csv(opt['mageck'], sep='\t')[['id', 'neg|fdr']].rename(columns={'id': 'gene', 'neg|fdr': 'mageck_neg_fdr'})
        t['mageck_hit'] = t['mageck_neg_fdr'] < fdr
        frames.append(t)
    if opt['bagel']:
        t = pd.read_csv(opt['bagel'], sep='\t')[['GENE', 'BF']].rename(columns={'GENE': 'gene', 'BF': 'bagel_bf'})
        t['bagel_hit'] = t['bagel_bf'] > bf
        frames.append(t)
    if opt['drugz']:
        t = pd.read_csv(opt['drugz'], sep='\t')[['GENE', 'fdr_synth']].rename(columns={'GENE': 'gene', 'fdr_synth': 'drugz_synth_fdr'})
        t['drugz_hit'] = t['drugz_synth_fdr'] < fdr
        frames.append(t)
    if len(frames) < 2:
        sys.exit('give at least two of mageck=, bagel=, drugz=: a consensus of one method is not a consensus')
    merged = frames[0]
    for t in frames[1:]:
        merged = merged.merge(t, on='gene', how='outer')
    flags = [c for c in merged.columns if c.endswith('_hit')]
    merged[flags] = merged[flags].fillna(False).astype(bool)
    merged['n_hits'] = merged[flags].sum(axis=1)
    merged['tier'] = merged['n_hits'].map({3: 'Tier-1', 2: 'Tier-2', 1: 'Tier-3'}).fillna('none')
    merged.sort_values(['n_hits', 'gene'], ascending=[False, True]).to_csv(out, index=False)
    print(merged['tier'].value_counts().to_string(), f'| methods: {len(flags)} | wrote {out}', sep='\n')


if __name__ == '__main__':
    main(sys.argv[1:])
