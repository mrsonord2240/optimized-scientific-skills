#!/usr/bin/env python3
'''PTM-SEA helpers for ssGSEA2.0: build the site-level GCT PTMsigDB expects, and read the scores back.

PTM-SEA is ssGSEA run on site identifiers. Each row is ONE localized site, named by its +/-7 flanking sequence
plus '-p' (e.g. AALRQLRSPRRAQAP-p): the central 15 residues of MaxQuant's 31-residue 'Sequence window'.

write-gct  table with columns `seq_window` (MaxQuant 'Sequence window', 31 aa) and `stat` (a SIGNED per-site
           statistic: log2FC or the moderated t from the protein-adjusted model, not intensities) -> GCT file.
           Rows whose window is not 31 aa are dropped; duplicate ids are averaged (a choice: report it).
read       ssGSEA2.0 output prefix -> table of NES and FDR per signature, sorted by NES.

The ssgsea-cli.R call itself is in references/ptm-sea.md.
Usage: python scripts/ptmsea.py write-gct --sites sites.tsv --out sites.gct
       python scripts/ptmsea.py read --prefix ptmsea/run --out ptmsea_scores.csv
Import: `from ptmsea import write_ptmsea_gct, read_ptmsea`.
'''
import argparse

import pandas as pd


def write_ptmsea_gct(sites, path):
    # sites: one row per LOCALIZED site, with `seq_window` (MaxQuant "Sequence window": 31 aa, 15 each side)
    # and `stat` (signed statistic). PTMsigDB wants the central 15 (+/-7) plus '-p'.
    t = sites[sites['seq_window'].str.len() == 31].copy()
    t['id'] = t['seq_window'].str[8:23] + '-p'
    t = t.groupby('id', as_index=False)['stat'].mean()   # ids must be unique; averaging duplicates is a choice, report it
    with open(path, 'w', newline='\n') as f:
        f.write(f'#1.2\n{len(t)}\t1\nName\tDescription\tsample1\n')
        for i, v in zip(t['id'], t['stat']):
            f.write(f'{i}\t{i}\t{v:.6f}\n')
    return t


def read_ptmsea(prefix='ptmsea/run'):
    gct = lambda kind: pd.read_csv(f'{prefix}-{kind}.gct', sep='\t', skiprows=2, index_col=0).iloc[:, -1]
    res = pd.DataFrame({'NES': gct('scores'), 'FDR': gct('fdr-pvalues')})
    return res.sort_values('NES', ascending=False)   # the run also writes -pvalues.gct and -combined.gct (signature size and overlap)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    w = sub.add_parser('write-gct')
    w.add_argument('--sites', required=True, help='TSV/CSV with seq_window and stat columns')
    w.add_argument('--out', default='sites.gct')
    r = sub.add_parser('read')
    r.add_argument('--prefix', default='ptmsea/run')
    r.add_argument('--out')
    a = ap.parse_args()
    if a.cmd == 'write-gct':
        sites = pd.read_csv(a.sites, sep=None, engine='python')
        t = write_ptmsea_gct(sites, a.out)
        print(f'{len(t)} unique site ids written to {a.out} ({len(sites) - len(t)} rows dropped or merged)')
    else:
        res = read_ptmsea(a.prefix)
        if a.out:
            res.to_csv(a.out)
        print(res.head(10).to_string())


if __name__ == '__main__':
    main()
