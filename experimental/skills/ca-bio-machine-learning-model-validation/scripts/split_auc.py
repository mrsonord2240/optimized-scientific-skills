#!/usr/bin/env python3
'''ROC AUC of a fixed pipeline under leave-one-site-out or forward-in-time splits.

Pipeline (selection, scaling, logistic regression) is refit on each training split. Use --site to test
transport to a new site, --time to test on later samples.

Inputs:  CSV with one row per sample, a binary label column, a site or time column, optional ID columns,
         and numeric feature columns (everything else).
Usage:   python split_auc.py data.csv --label label --site site [--drop sample_id] [--k 50] [--C 1.0]
         python split_auc.py data.csv --label label --time collection_date [--splits 5]
Output:  one AUC per held-out site (or per forward split), their mean and range.
'''
import argparse

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut, TimeSeriesSplit

from cv_auc import default_pipeline

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('csv')
    parser.add_argument('--label', required=True)
    parser.add_argument('--site', help='column naming the site; each site is held out once')
    parser.add_argument('--time', help='column giving time order (numeric or parseable date)')
    parser.add_argument('--drop', nargs='*', default=[], help='non-feature columns to ignore, e.g. sample_id')
    parser.add_argument('--k', type=int, default=50)
    parser.add_argument('--C', type=float, default=1.0)
    parser.add_argument('--splits', type=int, default=5, help='forward splits for --time')
    args = parser.parse_args()
    if bool(args.site) == bool(args.time):
        parser.error('give exactly one of --site and --time')

    data = pd.read_csv(args.csv)
    split_col = args.site or args.time
    if args.time:
        try:
            order = pd.to_datetime(data[split_col])
        except (ValueError, TypeError):
            order = data[split_col]
        data = data.assign(_order=order).sort_values('_order', kind='stable').drop(columns='_order')
    non_features = [args.label, split_col, *args.drop]
    features = data.drop(columns=non_features).select_dtypes('number')
    X, y = features.to_numpy(float), data[args.label].to_numpy()
    pipe = default_pipeline(min(args.k, features.shape[1]), args.C)

    if args.site:
        splits = LeaveOneGroupOut().split(X, y, data[args.site].to_numpy())
        names = lambda test: str(data[args.site].to_numpy()[test][0])
    else:
        splits = TimeSeriesSplit(n_splits=args.splits).split(X)
        names = lambda test: f'rows {test[0]}-{test[-1]}'

    aucs = []
    print(f'{args.csv}: {len(y)} rows, {features.shape[1]} features')
    for train, test in splits:
        if len(set(y[test])) < 2 or len(set(y[train])) < 2:
            print(f'{names(test)}: skipped (one class in train or test)')
            continue
        score = clone(pipe).fit(X[train], y[train]).decision_function(X[test])
        aucs.append(roc_auc_score(y[test], score))
        print(f'{names(test)}: n_test={len(test)} AUC={aucs[-1]:.2f}')
    if aucs:
        print(f'Mean AUC = {np.mean(aucs):.2f} (range {min(aucs):.2f}-{max(aucs):.2f} over {len(aucs)} '
              f'{"sites" if args.site else "forward splits"}; selection and scaling refit in every training split)')
