#!/usr/bin/env python
"""Dose-consistent drugZ hits.

Purpose: after running drugZ once per dose against the same vehicle, keep genes whose normZ has
the same sign at every dose and that pass FDR at the top dose. Reports whether |normZ| grows
with dose (supporting evidence only).
Inputs:  one drugZ output file per dose (tab-separated, columns GENE, normZ, fdr_synth, fdr_supp),
         given as DOSE=FILE pairs, and the name of the highest dose.
Usage:   python scripts/dose_consistent_hits.py --top-dose high low=drugz_low.txt mid=drugz_mid.txt high=drugz_high.txt [--fdr 0.05] [--direction synth|supp] [--out hits.tsv]
Import:  from dose_consistent_hits import dose_consistent_hits
"""
import argparse

import pandas as pd


def dose_consistent_hits(dose_files, top_dose, fdr=0.05, direction='synth'):
    """dose_files: {'low': 'drugz_low.txt', 'mid': ..., 'high': ...}; top_dose: key of the highest dose.

    A hit is dose-consistent if normZ has the same sign at every dose and passes FDR at the top dose.
    """
    frames = {d: pd.read_csv(f, sep='\t').set_index('GENE') for d, f in dose_files.items()}
    normz = pd.DataFrame({d: f['normZ'] for d, f in frames.items()}).dropna()
    sign_ok = (normz.gt(0).all(axis=1)) | (normz.lt(0).all(axis=1))
    top = frames[top_dose]
    passes = top['fdr_%s' % direction] < fdr
    hits = normz[sign_ok & passes.reindex(normz.index).fillna(False)].copy()
    hits['normZ_top_dose'] = top['normZ'].reindex(hits.index)
    hits['monotonic'] = (normz.abs().diff(axis=1).iloc[:, 1:] >= 0).all(axis=1)  # |normZ| grows with dose
    return hits.sort_values('normZ_top_dose')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('doses', nargs='+', metavar='DOSE=FILE', help='dose label and its drugZ output, in low-to-high order')
    ap.add_argument('--top-dose', required=True, help='label of the highest dose')
    ap.add_argument('--fdr', type=float, default=0.05)
    ap.add_argument('--direction', choices=['synth', 'supp'], default='synth')
    ap.add_argument('--out', help='write the hit table here (tab-separated); default: print')
    a = ap.parse_args()
    files = dict(x.split('=', 1) for x in a.doses)
    if a.top_dose not in files:
        ap.error('--top-dose %r is not one of the DOSE=FILE labels %s' % (a.top_dose, list(files)))
    hits = dose_consistent_hits(files, a.top_dose, a.fdr, a.direction)
    if a.out:
        hits.to_csv(a.out, sep='\t')
    else:
        print(hits.to_string())
    print('%d dose-consistent %s hits' % (len(hits), a.direction))
