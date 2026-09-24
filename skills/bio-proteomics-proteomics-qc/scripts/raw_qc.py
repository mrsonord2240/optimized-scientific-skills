#!/usr/bin/env python
# Purpose: inspect RAW (un-normalized) per-sample signal and strip contaminant rows BEFORE normalizing.
# Inputs:  a MaxQuant proteinGroups.txt, a sample sheet CSV (sample + group columns).
#          Un-normalized columns: MaxQuant "Intensity <sample>" (NOT "LFQ intensity <sample>").
# Usage:   python scripts/raw_qc.py proteinGroups.txt sample_annotation.csv [--group-col condition]
#          [--prefix "Intensity "]
#          or: sys.path.insert(0, "scripts"); from raw_qc import raw_sample_qc, strip_contaminant_rows
# Output:  per-sample n_quantified, total_signal, fold_total_vs_group, ids_vs_group, flag, loading_rule.
import argparse

import numpy as np
import pandas as pd

contaminant_flags = ['Potential contaminant', 'Reverse', 'Only identified by site']

def strip_contaminant_rows(protein_groups):
    keep = pd.Series(True, index=protein_groups.index)
    for col in contaminant_flags:
        match = next((c for c in protein_groups.columns if c.lower() == col.lower()), None)  # MaxQuant casing varies by version -- match case-insensitively
        if match is not None:
            keep &= protein_groups[match].fillna('') != '+'  # MaxQuant marks flagged rows with a literal '+'
    return protein_groups[keep]

def raw_sample_qc(raw_intensities, sample_groups):
    raw = raw_intensities.replace(0, np.nan)  # MaxQuant/DIA-NN write missing as 0; zeros are NOT quantified
    qc = pd.DataFrame({
        'n_quantified': raw.notna().sum(),
        'total_signal': raw.sum(),
        'median_intensity': raw.median(),
        'missing_pct': 100 * raw.isna().sum() / len(raw)})
    group = sample_groups.reindex(qc.index)
    sizes = group.value_counts()
    singletons = sorted(sizes[sizes < 2].index.astype(str))
    if singletons:
        # A group of one IS its own median: fold and ids would be exactly 1.000 and the flag
        # could never fire. Say so and fall back to the all-sample baseline for those samples.
        print(f'WARNING: single-sample group(s) {singletons}: the within-group loading rule cannot '
              'fire there, so those samples are compared to the ALL-sample median instead. A loading '
              'difference that tracks condition is NOT detectable in a design with no replicates.')
        baseline = group.where(group.map(sizes) >= 2, 'ALL')
    else:
        baseline = group
    qc['baseline'] = baseline
    qc['loading_rule'] = np.where(baseline == 'ALL', 'fallback_all_samples', 'within_group')
    qc['fold_total_vs_group'] = qc['total_signal'] / qc.groupby(baseline)['total_signal'].transform('median')
    qc['ids_vs_group'] = qc['n_quantified'] / qc.groupby(baseline)['n_quantified'].transform('median')
    qc['flag'] = (qc['fold_total_vs_group'] <= 0.5) | (qc['ids_vs_group'] < 0.8)  # >=2x low total, or >20% fewer IDs
    return qc

def contaminant_fraction(protein_groups, intensity_cols, flag_col='Potential contaminant'):
    flagged = protein_groups[flag_col].fillna('') == '+'
    raw = protein_groups[intensity_cols].replace(0, np.nan)
    return 100 * raw[flagged].sum() / raw.sum()  # percent of summed raw intensity, per sample


def main():
    ap = argparse.ArgumentParser(description=__doc__ or 'Raw per-sample QC before normalizing')
    ap.add_argument('protein_groups')
    ap.add_argument('sample_sheet')
    ap.add_argument('--sample-col', default='sample')
    ap.add_argument('--group-col', default='condition')
    ap.add_argument('--prefix', default='Intensity ', help='raw intensity column prefix (note the trailing space)')
    a = ap.parse_args()
    pg = pd.read_csv(a.protein_groups, sep='\t', low_memory=False)
    sheet = pd.read_csv(a.sample_sheet).set_index(a.sample_col)
    pg = strip_contaminant_rows(pg)
    cols = {c[len(a.prefix):]: c for c in pg.columns if c.startswith(a.prefix)}
    samples = [s for s in sheet.index if s in cols]
    raw = pg[[cols[s] for s in samples]]
    raw.columns = samples
    print(raw_sample_qc(raw, sheet[a.group_col]).round(3).to_string())


if __name__ == '__main__':
    main()
