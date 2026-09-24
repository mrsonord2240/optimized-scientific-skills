#!/usr/bin/env python
# Purpose: matrix-level QC: replicate correlation (log2), sample-swap check, CV on the linear scale,
#          missingness profile, completeness filter.
# Inputs:  a matrix (proteins x samples) and a sample sheet CSV; log2_intensities for correlation and
#          missingness, linear_intensities for CV. CLI reads MaxQuant proteinGroups.txt "LFQ intensity <sample>".
# Usage:   python scripts/matrix_metrics.py proteinGroups.txt sample_annotation.csv [--group-col condition]
#          or: sys.path.insert(0, "scripts"); from matrix_metrics import replicate_correlation, median_cv_linear
# Output:  DataFrames with a status column; filter status == 'measured' before summarizing.
import argparse
from itertools import combinations

import numpy as np
import pandas as pd

def replicate_correlation(log2_intensities, sample_groups):
    # DESIGN REQUIREMENT: >=2 samples in at least one group. Without it this silently returned
    # an empty table, which reads like "no reproducibility problem" when it means "not measured".
    sizes = sample_groups.value_counts()
    if (sizes >= 2).sum() == 0:
        raise ValueError(f'replicate_correlation needs >=2 samples in a group; sizes are '
                         f'{sizes.to_dict()}. With one run per condition there is no replicate '
                         'reproducibility to measure -- do not report "no outliers found".')
    corr = log2_intensities.corr(method='pearson')  # log2 first: Pearson on raw is a high-abundance artifact
    rows = []
    for group in sample_groups.unique():
        members = sample_groups[sample_groups == group].index
        if len(members) < 2:
            # a row per unchecked group, so the caller sees "not measured" in the data, not only on stdout
            print(f'WARNING: group {group!r} has no within-group pair; its sample is UNCHECKED here.')
            rows.append({'group': group, 's1': members[0], 's2': None, 'r': np.nan, 'status': 'not_measurable_n1'})
        for s1, s2 in combinations(members, 2):
            rows.append({'group': group, 's1': s1, 's2': s2, 'r': corr.loc[s1, s2], 'status': 'measured'})
    return pd.DataFrame(rows)  # filter status == 'measured' before summarizing r

def cross_group_correlation(log2_intensities, sample_groups, top_n=300):
    # Sample-swap check: does a sample match ANOTHER group better than its own? Within-group r alone
    # cannot show it (a swapped pair still correlates 0.93-0.96 with everything). Centre each protein
    # on its mean over samples, on the most variable proteins, so condition -- not shared abundance -- drives r.
    complete = log2_intensities.dropna(how='any')
    top = complete.loc[complete.var(axis=1).nlargest(top_n).index]
    corr = top.sub(top.mean(axis=1), axis=0).corr()
    rows = []
    for s in corr.columns:
        mean_r = {g: corr.loc[s, [x for x in corr.columns if sample_groups[x] == g and x != s]].mean()
                  for g in sample_groups.unique()}
        own = sample_groups[s]
        others = {g: r for g, r in mean_r.items() if pd.notna(r)}
        best = max(others, key=others.get)
        measurable = pd.notna(mean_r[own])
        rows.append({'sample': s, 'own_group': own, 'r_own': mean_r[own], 'best_group': best,
                     'r_best': others[best], 'status': 'measured' if measurable else 'not_measurable_n1',
                     'possible_swap': bool(best != own) if measurable else None})
    return pd.DataFrame(rows)

def median_cv_linear(linear_intensities, sample_groups):
    # DESIGN REQUIREMENT: >=2 samples per group. A one-sample group yields NaN, and a table of
    # NaNs is not "excellent precision" -- fail loudly instead of returning it.
    sizes = sample_groups.value_counts()
    if (sizes >= 2).sum() == 0:
        raise ValueError(f'median_cv_linear needs >=2 samples in a group; sizes are {sizes.to_dict()}. '
                         'CV is undefined with no replicates -- report "not measurable", not NaN.')
    rows = []
    for group in sample_groups.unique():
        members = sample_groups[sample_groups == group].index
        block = linear_intensities[members]
        per_protein_cv = block.std(axis=1) / block.mean(axis=1)  # base CV formula REQUIRES linear scale
        n1 = len(members) < 2
        if n1:
            print(f'WARNING: group {group!r} has {len(members)} sample(s); its CV is NaN (undefined), not low.')
        rows.append({'group': group, 'median_cv_pct': 100 * per_protein_cv.median(),
                     'status': 'not_measurable_n1' if n1 else 'measured'})
    return pd.DataFrame(rows)

def geometric_cv_from_log(log_intensities):
    sigma = log_intensities.std(axis=1) * np.log(2)  # convert log2 SD to natural-log SD
    return 100 * np.sqrt(np.expm1(sigma ** 2))  # gCV = sqrt(exp(sigma^2) - 1)

def missingness_profile(log2_intensities, n_bins=10):
    present_per_protein = log2_intensities.notna().mean(axis=1)
    mean_abundance = log2_intensities.mean(axis=1)
    abundance_bin = pd.qcut(mean_abundance, n_bins, duplicates='drop')
    # present fraction per mean-abundance bin: rising-with-abundance = MNAR, flat = MCAR
    return present_per_protein.groupby(abundance_bin, observed=True).mean()

def completeness_filter(log2_intensities, sample_groups, min_valid_frac=0.7):
    keep = pd.Series(False, index=log2_intensities.index)
    for group in sample_groups.unique():
        block = log2_intensities[sample_groups[sample_groups == group].index]
        keep |= block.notna().mean(axis=1) >= min_valid_frac  # valid in >=70% of >=1 condition
    return log2_intensities[keep]


def main():
    ap = argparse.ArgumentParser(description='Matrix-level QC metrics')
    ap.add_argument('protein_groups')
    ap.add_argument('sample_sheet')
    ap.add_argument('--sample-col', default='sample')
    ap.add_argument('--group-col', default='condition')
    ap.add_argument('--prefix', default='LFQ intensity ')
    a = ap.parse_args()
    pg = pd.read_csv(a.protein_groups, sep='\t', low_memory=False)
    sheet = pd.read_csv(a.sample_sheet).set_index(a.sample_col)
    groups = sheet[a.group_col]
    linear = pd.DataFrame({s: pg[a.prefix + s] for s in groups.index}).replace(0, np.nan)
    log2 = np.log2(linear)
    print(replicate_correlation(log2, groups).round(3).to_string())
    print(cross_group_correlation(log2, groups).round(3).to_string())
    print(median_cv_linear(linear, groups).round(2).to_string())
    print(missingness_profile(log2).round(3).to_string())
    print('complete rows after filter:', len(completeness_filter(log2, groups)), 'of', len(log2))


if __name__ == '__main__':
    main()
