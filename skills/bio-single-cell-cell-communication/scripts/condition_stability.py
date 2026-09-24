#!/usr/bin/env python
"""Null-split stability check for a two-condition LIANA comparison.

Purpose: each condition has fewer cells than the pooled data, so sampling noise alone flags
ligand-receptor pairs as gained/lost. This script computes the real gained/lost set, then permutes
the condition labels WITHIN each cell type (same cell numbers and composition), repeats the
comparison n-null times, and reports the real count against the null range. Pairs that are also
gained/lost in >=50% of null splits are flagged noise-prone.

Inputs: AnnData (.h5ad) with log-normalized .X, a cell-type column and a two-level condition column.
Output: printed summary; with --out, a TSV of gained/lost pairs with their null frequency.
Usage:  python condition_stability.py adata.h5ad --cond-a control --cond-b stimulated \
            [--groupby cell_type] [--condition condition] [--n-null 10] [--n-perms 1000] \
            [--seed 1337] [--out stability.tsv]
"""
import argparse
from collections import Counter

import liana as li
import numpy as np
import scanpy as sc


def robust_set(ad, groupby, n_perms):
    li.mt.rank_aggregate(ad, groupby=groupby, resource_name='consensus',
                         expr_prop=0.1, use_raw=False, n_perms=n_perms, verbose=False)
    res = ad.uns['liana_res']
    sig = res[(res['specificity_rank'] < 0.05) & (res['magnitude_rank'] < 0.05)]
    return set(zip(sig['source'], sig['target'], sig['ligand_complex'], sig['receptor_complex']))


def gained_lost(ad, groupby, condition, cond_a, cond_b, n_perms):
    a, b = (robust_set(ad[ad.obs[condition] == c].copy(), groupby, n_perms) for c in (cond_a, cond_b))
    return b - a, a - b   # gained in cond_b, lost from cond_a


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('h5ad')
    p.add_argument('--cond-a', required=True, help='reference condition')
    p.add_argument('--cond-b', required=True, help='comparison condition')
    p.add_argument('--groupby', default='cell_type')
    p.add_argument('--condition', default='condition')
    p.add_argument('--n-null', type=int, default=10)
    p.add_argument('--n-perms', type=int, default=1000)
    p.add_argument('--seed', type=int, default=1337)
    p.add_argument('--out')
    args = p.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    gained, lost = gained_lost(adata, args.groupby, args.condition, args.cond_a, args.cond_b, args.n_perms)

    rng = np.random.default_rng(args.seed)
    null_counts, pair_freq = [], Counter()
    for _ in range(args.n_null):
        null = adata.copy()
        null.obs[args.condition] = (null.obs.groupby(args.groupby, observed=True)[args.condition]
                                    .transform(lambda s: rng.permutation(s.values)).values)
        g, l = gained_lost(null, args.groupby, args.condition, args.cond_a, args.cond_b, args.n_perms)
        null_counts.append(len(g) + len(l))
        pair_freq.update(g | l)

    real = len(gained) + len(lost)
    print('real gained+lost:', real, '| null median:', np.median(null_counts),
          'range:', min(null_counts), max(null_counts))
    noise_prone = {q for q in gained | lost if pair_freq[q] / args.n_null >= 0.5}
    print('noise-prone pairs:', len(noise_prone), 'of', real)
    if min(null_counts) <= real <= max(null_counts):
        print('real count is inside the null range: not evidence of a condition effect')

    if args.out:
        with open(args.out, 'w', encoding='utf-8') as fh:
            fh.write('direction\tsource\ttarget\tligand\treceptor\tnull_freq\tnoise_prone\n')
            for tag, s in (('gained', gained), ('lost', lost)):
                for q in sorted(s):
                    fh.write('\t'.join([tag, *map(str, q), f'{pair_freq[q] / args.n_null:.2f}',
                                        str(q in noise_prone)]) + '\n')


if __name__ == '__main__':
    main()
