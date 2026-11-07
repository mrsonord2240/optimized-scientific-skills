import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Filter DIA-NN q-values and write log2 protein matrix')
parser.add_argument('report')
parser.add_argument('output')
parser.add_argument('--context', choices=['global', 'lib'], default='global')
args = parser.parse_args()
report = pd.read_parquet(args.report)
context = ['Global.Q.Value', 'Global.PG.Q.Value'] if args.context == 'global' else ['Lib.Q.Value', 'Lib.PG.Q.Value']
columns = ['Q.Value', 'PG.Q.Value'] + context
missing = set(columns + ['Protein.Group', 'Run', 'PG.MaxLFQ']) - set(report.columns)
if missing:
    parser.error('missing columns: ' + ', '.join(sorted(missing)))
filt = report.loc[(report[columns] <= 0.01).all(axis=1)]
matrix = filt.pivot_table(index='Protein.Group', columns='Run', values='PG.MaxLFQ', aggfunc='first')
np.log2(matrix.replace(0, np.nan)).to_csv(args.output, sep='\t')
print(f'input rows={len(report)}, passing rows={len(filt)}, matrix={matrix.shape}')
