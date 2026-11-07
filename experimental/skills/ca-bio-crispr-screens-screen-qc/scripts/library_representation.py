#!/usr/bin/env python
"""Per-sample library representation for a pooled CRISPR screen.

Purpose: detected sgRNAs, % zero, % low-count (<30 reads), median/p10/p90 and p90/p10 skew per sample,
         so the bottleneck can be localised between plasmid -> Day-0 -> endpoint.
Input:   MAGeCK count table (tab-separated; sgRNA index, Gene column, one numeric column per sample).
Usage:   python library_representation.py screen.count.txt [--out library_representation.tsv]
Import:  from library_representation import library_representation, stage_specific_thresholds
"""
import argparse
import pandas as pd
import numpy as np


def library_representation(counts_df):
    '''Per-sample library coverage diagnostics.
    counts_df: rows = sgRNAs, columns = samples (numeric counts).'''
    out = pd.DataFrame(index=counts_df.columns)
    out['n_sgrnas_detected'] = (counts_df > 0).sum()
    out['pct_zero'] = (counts_df == 0).sum() / len(counts_df) * 100
    out['pct_lowcount'] = (counts_df < 30).sum() / len(counts_df) * 100
    out['median_count'] = counts_df.median()
    out['p10_count'] = counts_df.quantile(0.10)
    out['p90_count'] = counts_df.quantile(0.90)
    out['skew_ratio'] = out['p90_count'] / out['p10_count'].replace(0, np.nan)
    return out


def stage_specific_thresholds():
    '''Stage conventions: Joung 2017 (zero-count, skew) + MAGeCK-VISPR (Gini).'''
    return {
        'plasmid':  {'pct_zero_max': 0.5, 'skew_max': 2.0, 'gini_max': 0.10},   # skew 2.0 is a stricter modern convention; Joung 2017 states <10
        'day_0':    {'pct_zero_max': 1.0, 'skew_max': 2.5, 'gini_max': 0.12},
        'endpoint': {'pct_zero_max': 5.0, 'skew_max': 10.0, 'gini_max': 0.30},
    }


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('counts', help='MAGeCK count table (screen.count.txt)')
    ap.add_argument('--out', help='write the table to this TSV as well')
    args = ap.parse_args()
    counts = pd.read_csv(args.counts, sep='\t', index_col=0)
    if 'Gene' in counts.columns:
        counts = counts.drop(columns='Gene')
    table = library_representation(counts)
    print(table.round(3).to_string())
    if args.out:
        table.to_csv(args.out, sep='\t')
