#!/usr/bin/env python3
# Reference: BAGEL2 2.0 (hart-lab/bagel, build 115), pandas 2.2+
#
# Purpose: classify genes from BAGEL2 `bf` output as essential / neutral (/ tumor_suppressor).
#   Tumor-suppressor calls are made only for enrichment-type screens (--screen-type enrichment|both);
#   the default 'dropout' never calls them. Assay-control pseudo-genes are dropped first.
# Inputs:  bayes_factor.txt (tab-separated, columns GENE and BF, from `BAGEL.py bf`).
# Output:  table sorted by BF descending with a `call` column (stdout, or -o FILE).
# Usage:   python scripts/interpret_bagel.py bayes_factor.txt [--screen-type dropout|enrichment|both]
#              [--bf-essential 6] [--bf-tumor-suppressor -6] [--controls LacZ,luciferase,EGFP] [-o calls.tsv]
#          or: from interpret_bagel import interpret_bagel

import argparse
import sys
import warnings

import pandas as pd

# Assay-control pseudo-genes spiked into CRISPR libraries (never real biology) --
# verified present and dominating the naive tumor-suppressor call on real HAP1 TKOv3
# output; extend this set to match your library's own controls.
ASSAY_CONTROLS = {'LacZ', 'luciferase', 'EGFP'}


def interpret_bagel(bf_path, bf_essential=6, bf_tumor_suppressor=-6,
                     screen_type='dropout', control_genes=ASSAY_CONTROLS,
                     tumor_suppressor_frac_warn=0.05):
    '''Classify genes from BAGEL2 BF output.

    screen_type: 'dropout' (default) only calls `essential`. Tumor-suppressor calls
    require screen_type='enrichment' or 'both' -- per the Failure Modes section of SKILL.md,
    a pure dropout screen's negative-BF genes are noise, not tumor suppressors.
    '''
    df = pd.read_csv(bf_path, sep='\t')
    df = df[~df['GENE'].isin(control_genes)].copy()   # drop assay-control pseudo-genes
    df['call'] = 'neutral'
    df.loc[df['BF'] > bf_essential, 'call'] = 'essential'
    if screen_type in ('enrichment', 'both'):
        df.loc[df['BF'] < bf_tumor_suppressor, 'call'] = 'tumor_suppressor'
        frac = (df['call'] == 'tumor_suppressor').mean()
        if frac > tumor_suppressor_frac_warn:
            warnings.warn(
                f"{frac:.1%} of genes flagged tumor_suppressor -- implausibly high; "
                "this usually means a dropout-only screen is being scored for "
                "enrichment. Re-check screen_type and BF<-6 calls against literature "
                "before reporting."
            )
    return df.sort_values('BF', ascending=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__ or 'Classify BAGEL2 bf output')
    ap.add_argument('bf_path', help='bayes_factor.txt from BAGEL.py bf')
    ap.add_argument('--screen-type', default='dropout', choices=['dropout', 'enrichment', 'both'])
    ap.add_argument('--bf-essential', type=float, default=6)
    ap.add_argument('--bf-tumor-suppressor', type=float, default=-6)
    ap.add_argument('--controls', default=','.join(sorted(ASSAY_CONTROLS)),
                    help='comma-separated assay-control pseudo-genes to drop')
    ap.add_argument('-o', '--out', help='write table here instead of stdout')
    a = ap.parse_args()
    df = interpret_bagel(a.bf_path, a.bf_essential, a.bf_tumor_suppressor, a.screen_type,
                         {g for g in a.controls.split(',') if g})
    if a.out:
        df.to_csv(a.out, sep='\t', index=False)
    else:
        df.to_csv(sys.stdout, sep='\t', index=False)
    print(df['call'].value_counts().to_string(), file=sys.stderr)


if __name__ == '__main__':
    main()
