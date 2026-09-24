"""Import a DIA-NN report.parquet as a filtered log2 protein-by-run matrix.

Usage: python load_diann.py report.parquet
"""
import sys
import numpy as np
import pandas as pd


REQUIRED = {'Q.Value', 'PG.Q.Value', 'Global.PG.Q.Value', 'Protein.Group', 'Run', 'PG.MaxLFQ'}


def load_diann(path):
    report = pd.read_parquet(path)
    missing = REQUIRED.difference(report.columns)
    if missing:
        raise ValueError(f"Not a supported DIA-NN report.parquet; missing {sorted(missing)}")
    report = report[(report['Q.Value'] <= 0.01) & (report['PG.Q.Value'] <= 0.01)
                    & (report['Global.PG.Q.Value'] <= 0.01)]
    matrix = report.pivot_table(index='Protein.Group', columns='Run', values='PG.MaxLFQ', aggfunc='first')
    return np.log2(matrix.replace(0, np.nan))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python load_diann.py report.parquet')
    matrix = load_diann(sys.argv[1])
    print(f'DIA-NN log2 matrix: {matrix.shape[0]} protein groups x {matrix.shape[1]} runs | -inf: {np.isinf(matrix.to_numpy()).any()}')
