#!/usr/bin/env python3
"""Guide-count QC for a pooled CRISPR screen: per-sample Gini (on ln(count+1), as `mageck count` reports it),
zero fraction, depth, skew, and within-condition replicate Pearson.

Usage: python qc.py <counts.txt> <out.tsv> [plasmid=<sample>] [min_depth=300] [pattern=<regex>]
  counts.txt   tab table: column 1 = guide id, column 2 = gene, then one RAW integer count column per sample
               (mageck count's `.count.txt`, not `.count_normalized.txt`)
  plasmid=     sample name(s), comma-separated: the cloned library pool, not Day-0 cells. Held to the plasmid Gini gate (<0.1); others use endpoint (<0.2)
  min_depth=   reads per sgRNA gate (default 300)
  pattern=     regex stripped from a sample name to get its condition (default removes a trailing _r1/_rep2/_A/_2, a letter after a digit (T18A), or an R<n> after a digit before an underscore (C902R1_P1D14)). No group of two or more: QC INCOMPLETE, exit 1
Output: out.tsv, one row per sample; replicate Pearson per condition printed. Exit 1 when a gate fails or replicates cannot be grouped.
Gates: Gini, reads per sgRNA, replicate Pearson. pct_guides_gt25 and skew_p90_p10 are reported, not gated.
Checked: pandas 3.0.5, numpy 2.5.3 (Python 3.12).
"""
import re
import sys
from itertools import combinations

import numpy as np
import pandas as pd

DEFAULT_PATTERN = r'(?:_(?:[Rr]ep|[Rr])?\d+$|_[A-Z]$|(?<=\d)[A-Z]$|(?<=\d)R\d+(?=_))'


def gini(x):
    # mageck count definition: all guides, on ln(count+1)
    x = np.sort(np.log(np.asarray(x, dtype=float) + 1))
    if x.size == 0 or x.sum() == 0:
        return np.nan
    i = np.arange(1, x.size + 1)
    return 2.0 * (i * x).sum() / x.size / x.sum() - 1 - 1.0 / x.size


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    path, out = argv[0], argv[1]
    opt = {'plasmid': '', 'min_depth': '300', 'pattern': DEFAULT_PATTERN}
    for item in argv[2:]:
        key, _, value = item.partition('=')
        if key not in opt:
            sys.exit(f'unknown option: {item}')
        opt[key] = value
    table = pd.read_csv(path, sep='\t', index_col=0)
    table = table.drop(columns=table.columns[0])  # gene column
    if table.shape[1] < 2:
        sys.exit('need at least two sample columns')
    if not np.allclose(table.values, np.round(table.values)):
        sys.exit('counts are not integers: use the raw .count.txt, not a normalized table')
    plasmid = [s for s in opt['plasmid'].split(',') if s]
    if not plasmid:
        print('WARNING: no plasmid= given: every sample is held to the endpoint Gini gate (<0.2), so a plasmid pool is not held to <0.1')
    missing = [s for s in plasmid if s not in table.columns]
    if missing:
        sys.exit(f'plasmid sample(s) not in the table: {missing}; columns are {list(table.columns)}')

    n = len(table)
    sample = pd.DataFrame({
        'pct_zero': (table == 0).sum() / n * 100,
        'gini': table.apply(gini),
        'reads_per_sgrna': table.sum() / n,
        'pct_guides_gt25': (table > 25).sum() / n * 100,
        'skew_p90_p10': table.apply(lambda c: np.percentile(c, 90) / max(np.percentile(c, 10), 1)),
    })
    sample['gini_gate'] = [0.1 if s in plasmid else 0.2 for s in sample.index]
    sample['pass'] = (sample['gini'] < sample['gini_gate']) & (sample['reads_per_sgrna'] >= float(opt['min_depth']))
    sample.to_csv(out, sep='\t')
    print(sample.round(3).to_string())

    logc = np.log10(table + 1).corr()
    groups = {}
    for col in table.columns:
        groups.setdefault(re.sub(opt['pattern'], '', col), []).append(col)
    failed = not sample['pass'].all()
    for cond, cols in groups.items():
        if len(cols) < 2:
            continue
        r = np.mean([logc.loc[a, b] for a, b in combinations(cols, 2)])
        print(f'replicate Pearson {cond} ({len(cols)} samples): {r:.3f}' + ('' if r >= 0.8 else '  BELOW 0.8'))
        failed |= r < 0.8
    ungrouped = not any(len(c) > 1 for c in groups.values())
    if ungrouped:
        print('replicate Pearson NOT CHECKED: no condition has two samples under this pattern; pass pattern=<regex> to group replicates')
    print('QC', 'FAIL' if failed else 'INCOMPLETE' if ungrouped else 'PASS', '| wrote', out)
    sys.exit(1 if failed or ungrouped else 0)


if __name__ == '__main__':
    main(sys.argv[1:])
