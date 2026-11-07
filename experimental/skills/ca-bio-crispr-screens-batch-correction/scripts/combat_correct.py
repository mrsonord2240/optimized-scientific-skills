#!/usr/bin/env python
"""ComBat batch correction for a CRISPR-screen count table, with a pre-fit guard.

Inputs : a tab-separated count table (first column = sgRNA id; extra non-sample columns such as
         Gene are carried through untouched) and a tab-separated metadata table (first column =
         sample name, plus a batch column and, optionally, a condition column).
Outputs: corrected counts (same rows and columns as the input) and the list of features that were
         left as raw counts because ComBat cannot fit them.
Usage  : python combat_correct.py counts.txt metadata.txt corrected.txt uncorrected.txt [--batch-col batch] [--condition-col condition]
Import : from combat_correct import combat_correct   # (corrected_df, uncorrected_index)
Needs  : pip install combat  (pyComBat, provides combat.pycombat), pandas, numpy
"""
import argparse
import numpy as np
import pandas as pd
from combat.pycombat import pycombat

def combat_correct(counts_df, batch_vector, condition_vector=None, verbose=True):
    '''ComBat on log2 counts with optional condition covariate (`mod`).
    Preserves condition signal while removing batch shifts.
    Returns (corrected_counts, uncorrected): the matrix (same rows and columns as the input) and
    the index of features that were left as raw counts.

    pycombat takes the matrix as a DataFrame (rows = features, columns = samples) and both
    `batch` and `mod` as plain lists of labels -- it one-hot-encodes `mod` itself, so passing a
    pre-encoded array fails inside pycombat.

    ComBat divides each feature by its pooled residual variance after batch (and `mod`) are
    regressed out. A feature with none left -- constant in every batch, or fully explained by
    batch + condition (e.g. counts 0,0,1,1 in each batch with condition 0,0,1,1) -- corrupts the
    shared empirical-Bayes prior, with no exception and exit code 0. Depending on the data the
    result is an all-NaN matrix (568,720 of 568,720 values on real TKOv3 counts) or the feature
    silently collapsing to a constant with no NaN at all. The pre-fit filter below is what
    prevents both; the NaN check after the fit is only a backstop. A feature that is constant in
    ONE batch but varies in another is safe and is corrected.
    '''
    data = pd.DataFrame(np.log2(counts_df.values + 1),
                        index=counts_df.index, columns=counts_df.columns)

    # Pre-fit filter: residual sum of squares after regressing out the design pycombat will use.
    design = [pd.get_dummies(pd.Series(list(batch_vector)), dtype=float)]
    if condition_vector is not None:
        design.append(pd.get_dummies(pd.Series(list(condition_vector)), drop_first=True, dtype=float))
    X = pd.concat(design, axis=1).to_numpy()
    resid = data.values - data.values @ (X @ np.linalg.pinv(X))
    usable = pd.Series((resid ** 2).sum(axis=1) > 1e-8, index=data.index)   # not == 0: round-off
    uncorrected = data.index[~usable]
    if verbose and len(uncorrected):
        print(f'ComBat: {len(uncorrected)} features have no variance left after batch and condition; '
              f'returned as raw counts (see the `uncorrected` return value)')

    if condition_vector is not None:
        corrected = pycombat(data[usable], list(batch_vector), mod=list(condition_vector))
    else:
        corrected = pycombat(data[usable], list(batch_vector))

    # Backstop only -- fail loudly rather than pass a dead matrix on.
    if corrected.isna().any().any():
        raise ValueError('ComBat returned NaN values despite the pre-fit filter; do not use the '
                         'output. Inspect the design (batch/condition labels) and the count matrix.')

    out = pd.DataFrame(np.power(2, corrected.values) - 1,
                       index=corrected.index, columns=corrected.columns).clip(lower=0)
    return out.reindex(counts_df.index).fillna(counts_df), uncorrected   # dropped features keep raw counts


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('counts')
    ap.add_argument('metadata')
    ap.add_argument('out_counts')
    ap.add_argument('out_uncorrected')
    ap.add_argument('--batch-col', default='batch')
    ap.add_argument('--condition-col', default='condition',
                    help='column passed as `mod`; use "" to run without a covariate')
    a = ap.parse_args()

    raw = pd.read_csv(a.counts, sep='	', index_col=0)
    meta = pd.read_csv(a.metadata, sep='	', index_col=0)
    samples = [c for c in raw.columns if c in meta.index]
    if not samples:
        raise SystemExit('No count-table column matches a metadata sample name')
    meta = meta.loc[samples]
    cond = list(meta[a.condition_col]) if a.condition_col else None
    corrected, uncorrected = combat_correct(raw[samples], list(meta[a.batch_col]), cond)
    extra = raw.drop(columns=samples)
    pd.concat([extra, corrected], axis=1).to_csv(a.out_counts, sep='	')
    pd.Series(uncorrected, name='uncorrected').to_csv(a.out_uncorrected, sep='	', index=False)
    print(f'{len(corrected)} features corrected or passed through; {len(uncorrected)} left uncorrected')


if __name__ == '__main__':
    main()
