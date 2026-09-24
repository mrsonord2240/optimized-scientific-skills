#!/usr/bin/env python
# Purpose: per-run Level-1 metrics from a DIA-NN report: RT vs Predicted.RT R^2, median FWHM, Quantity.Quality.
# Inputs:  DIA-NN report.parquet or report.tsv (columns Run, RT, Predicted.RT, FWHM, Quantity.Quality, Global.Q.Value).
# Usage:   python scripts/diann_level1.py report.parquet
#          or: sys.path.insert(0, "scripts"); from diann_level1 import diann_level1
# Output:  one row per run with rt_fit_r2, median_fwhm, fwhm_vs_median, flag.
import argparse

import numpy as np
import pandas as pd

def diann_level1(report):
    # report: DIA-NN report.tsv/parquet as a DataFrame (columns Run, RT, Predicted.RT, FWHM, Quantity.Quality, Global.Q.Value)
    rep = report[report['Global.Q.Value'] <= 0.01].dropna(subset=['RT', 'Predicted.RT', 'FWHM'])
    rows = []
    for run, d in rep.groupby('Run'):
        rows.append({'run': run, 'n_precursors': len(d),
                     'rt_fit_r2': np.corrcoef(d['RT'], d['Predicted.RT'])[0, 1] ** 2,
                     'median_fwhm': d['FWHM'].median(),
                     'median_quantity_quality': d['Quantity.Quality'].median()})  # reported only: no accepted cutoff
    out = pd.DataFrame(rows).set_index('run')
    out['fwhm_vs_median'] = out['median_fwhm'] / out['median_fwhm'].median()
    out['flag'] = (out['rt_fit_r2'] < 0.99) | (out['fwhm_vs_median'] > 1.25)
    return out


def main():
    ap = argparse.ArgumentParser(description='DIA-NN Level-1 run metrics')
    ap.add_argument('report')
    a = ap.parse_args()
    report = pd.read_parquet(a.report) if a.report.endswith('.parquet') else pd.read_csv(a.report, sep='\t')
    print(diann_level1(report).round(4).to_string())


if __name__ == '__main__':
    main()
