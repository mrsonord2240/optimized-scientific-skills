#!/usr/bin/env python3
"""Two-condition screen: run `mageck test` (RRA), write hit tables and a volcano plot.

Usage: python rra.py <counts.txt> <prefix> treatment=<s1,s2> control=<s3> [norm=median] [fdr=0.05] [lfc=0.5]
  counts.txt   mageck count table (raw counts). Cancer line: the CRISPRcleanR-corrected table.
  treatment=   sample labels of the endpoint arm (required; there is no default)
  control=     sample labels of the baseline: plasmid, Day0, or vehicle for a drug screen (required)
  norm=        median (default) or control (with ctrl_genes=<file> of non-targeting guide genes)
  fdr=, lfc=   hit cutoffs; a hit needs fdr below `fdr` and |lfc| of at least `lfc`. lfc=0 drops no row.
Output: <prefix>.gene_summary.txt, <prefix>_negative_hits.csv, <prefix>_positive_hits.csv, <prefix>_volcano.png
Checked: MAGeCK 0.5.9.5, pandas 3.0.5, matplotlib 3.11.1.
"""
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    path, prefix = argv[0], argv[1]
    opt = {'treatment': '', 'control': '', 'norm': 'median', 'fdr': '0.05', 'lfc': '0.5', 'ctrl_genes': ''}
    for item in argv[2:]:
        key, _, value = item.partition('=')
        if key not in opt:
            sys.exit(f'unknown option: {item}')
        opt[key] = value
    if not opt['treatment'] or not opt['control']:
        sys.exit('treatment= and control= are required: the baseline decides every LFC '
                 '(plasmid/Day0 for dropout, vehicle for a drug screen)')
    header = pd.read_csv(path, sep='\t', nrows=0).columns
    treatment, control = opt['treatment'].split(','), opt['control'].split(',')
    for label in treatment + control:
        if label not in header:
            sys.exit(f'sample "{label}" is not a column of {path}: {list(header[2:])}')
    if set(treatment) & set(control):
        sys.exit('a sample is in both treatment and control')
    if opt['norm'] == 'control' and not opt['ctrl_genes']:
        sys.exit('norm=control needs ctrl_genes=<file listing the non-targeting control genes>')
    exe = shutil.which('mageck')
    launcher = [exe]
    if exe is None:  # Windows installs `mageck` as an extensionless Python script
        for folder in os.environ.get('PATH', '').split(os.pathsep):
            if os.path.isfile(os.path.join(folder, 'mageck')):
                exe, launcher = os.path.join(folder, 'mageck'), [sys.executable, os.path.join(folder, 'mageck')]
                break
    if exe is None:
        sys.exit('mageck is not on PATH (see references/install.md)')

    cmd = launcher + ['test', '--count-table', path, '--treatment-id', opt['treatment'], '--control-id', opt['control'],
           '--norm-method', opt['norm'], '--output-prefix', prefix]
    if opt['norm'] == 'control':
        cmd += ['--control-gene', opt['ctrl_genes']]
    subprocess.run(cmd, check=True, capture_output=True)

    gs = pd.read_csv(f'{prefix}.gene_summary.txt', sep='\t')
    fdr, lfc = float(opt['fdr']), float(opt['lfc'])
    neg = gs[(gs['neg|fdr'] < fdr) & (gs['neg|lfc'].abs() >= lfc) & (gs['neg|lfc'] < 0)].sort_values('neg|rank')
    pos = gs[(gs['pos|fdr'] < fdr) & (gs['pos|lfc'].abs() >= lfc) & (gs['pos|lfc'] > 0)].sort_values('pos|rank')
    neg.to_csv(f'{prefix}_negative_hits.csv', index=False)
    pos.to_csv(f'{prefix}_positive_hits.csv', index=False)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 6))
    for side, colour in (('neg', 'tab:red'), ('pos', 'tab:blue')):
        y = -np.log10(gs[f'{side}|fdr'].clip(lower=1e-10))
        sig = gs[f'{side}|fdr'] < fdr
        ax.scatter(gs.loc[~sig, f'{side}|lfc'], y[~sig], c='lightgray', s=8, alpha=0.5)
        ax.scatter(gs.loc[sig, f'{side}|lfc'], y[sig], c=colour, s=14, alpha=0.7)
    ax.axhline(-np.log10(fdr), ls='--', c='black', lw=0.5)
    ax.set_xlabel('log2 fold change')
    ax.set_ylabel('-log10(FDR)')
    fig.savefig(f'{prefix}_volcano.png', dpi=150)
    print(f'genes: {len(gs)} | negative hits: {len(neg)} | positive hits: {len(pos)} | wrote {prefix}.gene_summary.txt')


if __name__ == '__main__':
    main(sys.argv[1:])
