import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('output')
parser.add_argument('--heavy', default='Intensity H')
parser.add_argument('--light', default='Intensity L')
args = parser.parse_args()
table = pd.read_csv(args.input, sep='\t')
heavy = np.asarray(table[args.heavy], dtype=float)
light = np.asarray(table[args.light], dtype=float)
h, l = heavy > 0, light > 0
with np.errstate(divide='ignore', invalid='ignore'):
    table['log2_HL'] = np.where(h & l, np.log2(heavy / light), np.nan)
table['presence'] = np.select([h & l, h, l], ['both', 'H-only', 'L-only'], default='none')
table.to_csv(args.output, sep='\t', index=False)
