'''Expand MaxQuant Phospho (STY)Sites multiplicity and filter to class I localization.

Self-contained: writes a tiny synthetic site table to a tempdir, processes it, cleans up.
The point is the multiplicity expansion (Intensity <run>___1/___2/___3, THREE underscores)
and the protein-adjustment caveat - quantifying on the collapsed base Intensity mixes
phospho-states and can mask a switch between singly- and doubly-phosphorylated forms.
The table carries both the aggregated Intensity___n columns (sum over runs) and the
per-run columns, as MaxQuant writes them; only the per-run columns are melted.
'''
# Reference: pandas 2.2+, numpy 1.26+ | Verify API if version differs
import os
import re
import tempfile
import pandas as pd
import numpy as np

CLASS_I_PROB = 0.75   # Olsen 2006 class-I convention; comparability standard, not a calibrated FLR
RUNS = ['C1', 'T1']

synthetic = pd.DataFrame({
    'Gene names': ['AKT1;AKT1b', 'MAPK1', np.nan, 'GSK3B'],
    'Protein': ['P31749', 'P28482', 'Q9Y6K9', 'P49841'],
    'Amino acid': ['S', 'T', 'Y', 'S'],
    'Position': [473, 185, 99, 9],
    'Localization prob': [0.99, 0.75, 0.61, 0.95],   # 0.75 is class I (>=); 0.61 is class II, filtered out
    'Sequence window': ['_' * 31, '_' * 31, '_' * 31, '_' * 31],
    'Reverse': ['', '', '', ''],
    'Potential contaminant': ['', '', '', ''],
    # per-run columns; GSK3B S9 switches from singly (C1) to doubly phosphorylated (T1)
    'Intensity C1___1': [6.0e7, 3.0e7, 1.0e7, 8.0e7],
    'Intensity C1___2': [1.5e7, 0.0, 0.0, 1.0e7],
    'Intensity T1___1': [4.0e7, 2.0e7, 1.0e7, 2.0e7],
    'Intensity T1___2': [1.5e7, 0.0, 0.0, 7.0e7],
})
for n in ['1', '2']:   # aggregated sum-over-runs columns, as MaxQuant also writes them
    synthetic[f'Intensity___{n}'] = synthetic[[f'Intensity {r}___{n}' for r in RUNS]].sum(axis=1)

with tempfile.TemporaryDirectory() as tmp:
    synthetic.to_csv(os.path.join(tmp, 'Phospho (STY)Sites.txt'), sep='\t', index=False)

    # current MaxQuant writes a space in the modification name; accept either form
    path = next(p for p in (os.path.join(tmp, f) for f in ['Phospho (STY)Sites.txt', 'Phospho(STY)Sites.txt']) if os.path.exists(p))
    phospho = pd.read_csv(path, sep='\t', low_memory=False)
    print(f'Total sites: {len(phospho)}')

    contaminant_col = 'Potential contaminant' if 'Potential contaminant' in phospho.columns else 'Contaminant'
    phospho = phospho[(phospho['Reverse'] != '+') & (phospho[contaminant_col] != '+')]

    # right=False so the bins match the >= filter: [0.75, 1] is class I
    phospho['loc_class'] = pd.cut(phospho['Localization prob'], bins=[0, 0.25, 0.5, 0.75, 1.0 + 1e-9],
                                  labels=['IV', 'III', 'II', 'I'], right=False)
    print('\nLocalization class distribution:')
    print(phospho['loc_class'].value_counts())

    confident = phospho[phospho['Localization prob'] >= CLASS_I_PROB].copy()
    print(f'\nClass I sites (prob >= {CLASS_I_PROB}): {len(confident)}')

    gene = confident['Gene names'].where(confident['Gene names'].notna(), confident['Protein'])
    confident['site_id'] = gene.str.split(';').str[0] + '_' + confident['Amino acid'] + confident['Position'].astype(int).astype(str)

    # per-run multiplicity columns only: the space after 'Intensity' excludes the aggregated Intensity___n
    mult_cols = [c for c in confident.columns if re.fullmatch(r'Intensity .+___[123]', c)]
    long = confident.melt(id_vars=['site_id', 'Amino acid', 'Position', 'Localization prob'], value_vars=mult_cols, var_name='run_multiplicity', value_name='intensity')
    long['multiplicity'] = long['run_multiplicity'].str.split('___').str[-1]
    long['run'] = long['run_multiplicity'].str.replace(r'___[123]$', '', regex=True).str.replace(r'^Intensity ', '', regex=True)
    long = long[long['intensity'] > 0].copy()
    long['log2_intensity'] = np.log2(long['intensity'])
    print(f'Runs in long table: {sorted(long["run"].unique())}')

    print('\nMultiplicity-resolved site observations (collapsing these would mask the GSK3B form switch):')
    print(long.pivot_table(index=['site_id', 'multiplicity'], columns='run', values='log2_intensity').round(2).to_string())

    print('\nNext step (not shown): feed the PTM site table AND a paired global proteome to')
    print('MSstatsPTM groupComparisonPTM and call only ADJUSTED.Model hits regulated.')
