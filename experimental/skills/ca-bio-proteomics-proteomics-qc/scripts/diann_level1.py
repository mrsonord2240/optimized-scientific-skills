#!/usr/bin/env python
# Purpose: per-run Level-1 metrics from a DIA-NN report: RT vs Predicted.RT R^2, median FWHM, Quantity.Quality.
# Inputs:  DIA-NN report.parquet or report.tsv (columns Run, RT, Predicted.RT, FWHM, Quantity.Quality, Global.Q.Value).
# Usage:   python scripts/diann_level1.py report.parquet
#          or: sys.path.insert(0, "scripts"); from diann_level1 import diann_level1
# Output:  one row per run with rt_fit_r2, median_fwhm, fwhm_vs_median, flag.
import argparse
import sys

import numpy as np
import pandas as pd

def diann_level1(report, max_q=0.01, min_r2=0.99, max_fwhm_ratio=1.25):
    # report: DIA-NN report.tsv/parquet as a DataFrame (columns Run, RT, Predicted.RT, FWHM, Quantity.Quality, Global.Q.Value)
    rep = report[report['Global.Q.Value'] <= max_q].dropna(subset=['RT', 'Predicted.RT', 'FWHM'])
    rows = []
    for run, d in rep.groupby('Run'):
        rows.append({'run': run, 'n_precursors': len(d),
                     'rt_fit_r2': np.corrcoef(d['RT'], d['Predicted.RT'])[0, 1] ** 2,
                     'median_fwhm': d['FWHM'].median(),
                     'median_quantity_quality': d['Quantity.Quality'].median()})  # reported only: no accepted cutoff
    out = pd.DataFrame(rows).set_index('run')
    out['fwhm_vs_median'] = out['median_fwhm'] / out['median_fwhm'].median()
    out['flag'] = (out['rt_fit_r2'] < min_r2) | (out['fwhm_vs_median'] > max_fwhm_ratio)
    return out


def main():
    ap = argparse.ArgumentParser(description='DIA-NN Level-1 run metrics')
    ap.add_argument('report')
    ap.add_argument('--max-q', type=float, default=0.01)
    ap.add_argument('--min-r2', type=float, default=0.99, help='flag a run whose RT fit R^2 is below this')
    ap.add_argument('--max-fwhm-ratio', type=float, default=1.25, help='flag a run whose median FWHM exceeds this multiple of the across-run median')
    ap.add_argument('--out', default='diann_level1_results.tsv')
    a = ap.parse_args()
    report = pd.read_parquet(a.report) if a.report.endswith('.parquet') else pd.read_csv(a.report, sep='\t')
    need = ['Run', 'RT', 'Predicted.RT', 'FWHM', 'Quantity.Quality', 'Global.Q.Value']
    absent = [c for c in need if c not in report.columns]
    if absent:
        sys.exit(f'REFUSED: report lacks {absent}; this needs the main DIA-NN report (parquet or tsv), not a matrix file')
    out = diann_level1(report, a.max_q, a.min_r2, a.max_fwhm_ratio)
    if len(out) < 4:
        print(f'WARNING: {len(out)} runs; the across-run median is a weak baseline, trend over the queue before acting on a flag.')
    out.round(4).to_csv(a.out, sep='\t')
    print(out.round(4).to_string())
    print('wrote', a.out)


if __name__ == '__main__':
    main()
