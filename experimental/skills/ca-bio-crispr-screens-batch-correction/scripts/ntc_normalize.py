#!/usr/bin/env python3
import argparse
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('counts')
    parser.add_argument('controls')
    parser.add_argument('normalized')
    parser.add_argument('factors')
    args = parser.parse_args()
    counts = pd.read_csv(args.counts, sep='\t', index_col=0)
    with open(args.controls, encoding='utf-8') as handle:
        controls = {line.strip() for line in handle if line.strip()}
    selected = counts.index.isin(controls)
    if not selected.any():
        parser.error('No control IDs matched the count matrix')
    if selected.sum() < 500:
        print('Warning: fewer than 500 matched controls', file=sys.stderr)
    medians = counts.loc[selected].median(axis=0)
    if medians.isna().any() or (medians == 0).any():
        parser.error('Control median is missing or zero in at least one sample')
    factors = 1000 / medians
    (counts * factors).to_csv(args.normalized, sep='\t')
    factors.rename('scale_factor').to_csv(args.factors, sep='\t', header=True)

if __name__ == '__main__':
    main()
