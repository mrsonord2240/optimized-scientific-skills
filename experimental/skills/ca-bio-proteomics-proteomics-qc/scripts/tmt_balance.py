#!/usr/bin/env python
# Purpose: TMT channel balance within each plex, on RAW reporter intensities.
# Inputs:  one CSV per plex (rows = proteins, first column the protein id, one column per channel) of
#          un-normalized reporter intensities, i.e. before any channel or global normalization.
# Usage:   python scripts/tmt_balance.py plexA.csv plexB.csv [--max-fold 2] [--out tmt_balance_results.tsv]
#          Plex name = file stem. Never compare channels across plexes; use the reference channel for that.
# Output:  one TSV: plex, channel, total_signal, fold_vs_plex_median, investigate (beyond max-fold),
#          flag_severe (beyond 3.5x).
import argparse
import sys

import numpy as np
import pandas as pd


def tmt_channel_balance(plex_matrices, max_fold=2.0):
    # plex_matrices: {plex_name: DataFrame proteins x channels of RAW reporter intensities}
    rows = []
    for plex, m in plex_matrices.items():
        total = m.replace(0, np.nan).sum()
        for channel, fold in (total / total.median()).items():
            rows.append({'plex': plex, 'channel': channel, 'total_signal': total[channel],
                         'fold_vs_plex_median': fold})
    balance = pd.DataFrame(rows)
    balance['investigate'] = np.abs(np.log2(balance['fold_vs_plex_median'])) > np.log2(max_fold)
    balance['flag_severe'] = np.abs(np.log2(balance['fold_vs_plex_median'])) > np.log2(3.5)
    return balance


def main():
    ap = argparse.ArgumentParser(description='TMT channel balance per plex')
    ap.add_argument('plex_csv', nargs='+')
    ap.add_argument('--max-fold', type=float, default=2.0, help='investigate beyond this fold from the plex median')
    ap.add_argument('--out', default='tmt_balance_results.tsv')
    a = ap.parse_args()
    plexes = {}
    for f in a.plex_csv:
        m = pd.read_csv(f, index_col=0)
        if not all(np.issubdtype(t, np.number) for t in m.dtypes):
            sys.exit(f'REFUSED: {f} has non-numeric channel columns; pass the reporter-intensity table only')
        if m.min().min() < 0 or m.replace(0, np.nan).stack().median() < 100:
            sys.exit(f'REFUSED: {f} looks log-transformed or normalized (median reporter intensity < 100). '
                     'Balance is read on raw reporter intensities.')
        plexes[f.replace('\\', '/').split('/')[-1].rsplit('.', 1)[0]] = m
    out = tmt_channel_balance(plexes, a.max_fold)
    out.round(4).to_csv(a.out, sep='\t', index=False)
    print(out.round(3).to_string())
    print('wrote', a.out)


if __name__ == '__main__':
    main()
