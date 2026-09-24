#!/usr/bin/env python3
'''Kinase / writer motif enrichment around class-I phosphosites, against an experiment-matched background.

Foreground = the +/-7 windows (from MaxQuant 'Sequence window', a 31-mer centred on the site) of the class-I
sites, or of a chosen subset (--foreground-ids, e.g. the regulated-up sites). Background = every S/T (or Y)
window in the proteins identified in this experiment, central residue preserved. One-sided Fisher test per
(offset, residue), BH-adjusted across all tests. A whole-proteome or IUPAC-random background only rediscovers
the composition of disordered regions.

Inputs:
  --sites           MaxQuant 'Phospho (STY)Sites.txt' (the no-space 'Phospho(STY)Sites.txt' is tried too when a
                    directory is given)
  --fasta           search FASTA (UniProt headers; the accession is the second '|' field)
  --foreground-ids  optional text file, one site_id per line (Gene_S473, as in the multiplicity-expansion recipe);
                    default: every class-I site
  --residues        central residues of the background (default ST; use Y for tyrosine)
  --out             enrichment table (CSV)
Usage: python scripts/motif_enrichment.py --sites 'Phospho (STY)Sites.txt' --fasta search.fasta --out motif.csv
Import: `from motif_enrichment import matched_background, motif_enrichment` for the two helpers.
'''
import argparse
import os
from collections import Counter

import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control

CLASS_I_PROB = 0.75  # Olsen 2006 class-I convention; comparability standard, not a calibrated FLR
# 'Sequence window' is a 31-mer (+/-15) centered on the modified residue.
WINDOW_HALF = 7  # +/-7 flanking is the standard kinase-motif window


def read_fasta(path):
    fasta, acc = {}, None
    for line in open(path, encoding='utf-8'):
        line = line.strip()
        if line.startswith('>'):
            head = line[1:].split()[0]
            acc = head.split('|')[1] if head.count('|') >= 2 else head
            fasta[acc] = ''
        elif acc:
            fasta[acc] += line
    return fasta


def matched_background(fasta, accessions, residues='ST'):
    '''Every S/T (or Y) window in the IDENTIFIED proteins, central residue preserved.
    fasta: {accession: sequence}; accessions: proteins identified in this experiment.'''
    windows = []
    for acc in accessions:
        seq = fasta[acc]
        for i, aa in enumerate(seq):
            if aa in residues:
                windows.append(''.join(seq[j] if 0 <= j < len(seq) else '_' for j in range(i - WINDOW_HALF, i + WINDOW_HALF + 1)))
    return windows


def position_frequencies(windows):
    counts = {i: Counter() for i in range(-WINDOW_HALF, WINDOW_HALF + 1)}
    for w in windows:
        for offset, aa in zip(range(-WINDOW_HALF, WINDOW_HALF + 1), w):
            if aa not in '_X':
                counts[offset][aa] += 1
    return counts


def motif_enrichment(fg_windows, bg_windows):
    '''One-sided Fisher test per (position, residue), BH-adjusted across all tests.'''
    fg, bg = position_frequencies(fg_windows), position_frequencies(bg_windows)
    rows = []
    for offset in fg:
        if offset == 0:
            continue
        n_fg, n_bg = sum(fg[offset].values()), sum(bg[offset].values())
        for aa, k in fg[offset].items():
            _, p = fisher_exact([[k, n_fg - k], [bg[offset][aa], n_bg - bg[offset][aa]]], alternative='greater')
            rows.append({'offset': offset, 'aa': aa, 'fg': k, 'fg_total': n_fg, 'bg': bg[offset][aa], 'bg_total': n_bg, 'p': p})
    table = pd.DataFrame(rows)
    table['q'] = false_discovery_control(table['p'], method='bh')
    return table.sort_values('p')


def load_class_i_sites(path):
    '''Class-I sites of a MaxQuant Phospho (STY)Sites table, with Gene_<residue><position> site ids.'''
    if os.path.isdir(path):
        path = next(os.path.join(path, f) for f in ['Phospho (STY)Sites.txt', 'Phospho(STY)Sites.txt'] if os.path.exists(os.path.join(path, f)))
    phospho = pd.read_csv(path, sep='\t', low_memory=False)
    contaminant_col = 'Potential contaminant' if 'Potential contaminant' in phospho.columns else 'Contaminant'
    phospho = phospho[(phospho['Reverse'] != '+') & (phospho[contaminant_col] != '+')]
    phospho = phospho[phospho['Localization prob'] >= CLASS_I_PROB].copy()
    gene = phospho['Gene names'].where(phospho['Gene names'].notna(), phospho['Protein'])
    phospho['site_id'] = gene.str.split(';').str[0] + '_' + phospho['Amino acid'] + phospho['Position'].astype(int).astype(str)
    return phospho


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sites', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--foreground-ids')
    ap.add_argument('--residues', default='ST')
    ap.add_argument('--out', default='motif_enrichment.csv')
    a = ap.parse_args()

    phospho = load_class_i_sites(a.sites)
    phospho = phospho[phospho['Sequence window'].notna() & (phospho['Sequence window'].str.len() >= 31)]
    fg_rows = phospho
    if a.foreground_ids:
        ids = {l.strip() for l in open(a.foreground_ids, encoding='utf-8') if l.strip()}
        fg_rows = phospho[phospho['site_id'].isin(ids)]
    foreground = [w[15 - WINDOW_HALF: 16 + WINDOW_HALF] for w in fg_rows['Sequence window']]

    fasta = read_fasta(a.fasta)
    accessions = sorted({acc for p in phospho['Protein'] for acc in str(p).split(';') if acc in fasta})
    background = matched_background(fasta, accessions, residues=a.residues)
    print(f'foreground windows: {len(foreground)} | background windows: {len(background)} from {len(accessions)} proteins')

    table = motif_enrichment(foreground, background)
    table.to_csv(a.out, index=False)
    print(table.head(10).to_string(index=False))


if __name__ == '__main__':
    main()
