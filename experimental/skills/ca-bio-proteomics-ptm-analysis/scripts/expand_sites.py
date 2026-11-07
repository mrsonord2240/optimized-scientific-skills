import argparse
import re
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('sites')
parser.add_argument('output')
args = parser.parse_args()
phospho = pd.read_csv(args.sites, sep='\t', low_memory=False)
contaminant_col = 'Potential contaminant' if 'Potential contaminant' in phospho.columns else 'Contaminant'
phospho = phospho[(phospho['Reverse'] != '+') & (phospho[contaminant_col] != '+')]
confident = phospho[phospho['Localization prob'] >= 0.75].copy()
gene = confident['Gene names'].where(confident['Gene names'].notna(), confident['Protein'])
confident['site_id'] = gene.str.split(';').str[0] + '_' + confident['Amino acid'] + confident['Position'].astype(int).astype(str)
mult_cols = [c for c in confident.columns if re.fullmatch(r'Intensity .+___[123]', c)]
if not mult_cols:
    parser.error('no per-run multiplicity intensity columns found')
long = confident.melt(id_vars=['site_id', 'Amino acid', 'Position', 'Localization prob'], value_vars=mult_cols, var_name='run_multiplicity', value_name='intensity')
long['multiplicity'] = long['run_multiplicity'].str.split('___').str[-1]
long['run'] = long['run_multiplicity'].str.replace(r'___[123]$', '', regex=True).str.replace(r'^Intensity ', '', regex=True)
long = long[long['intensity'] > 0].copy()
long['log2_intensity'] = np.log2(long['intensity'])
long.to_csv(args.output, sep='\t', index=False)
