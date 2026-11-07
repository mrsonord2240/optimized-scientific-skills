#!/usr/bin/env python3
'''Nested cross-validated ROC AUC for a tuned pipeline (scaling, SelectKBest, logistic regression).

The inner loop picks the number of features k and the regularisation C; the outer loop grades the winning
configuration on folds it never saw. Every data-dependent step is refit inside each training fold.

Inputs:  CSV with one row per sample, a binary label column, optional ID and group columns, and numeric
         feature columns (everything else).
Usage:   python nested_cv.py data.csv --label label [--group patient_id] [--drop sample_id]
                             [--ks 10,50,200] [--Cs 0.01,0.1,1] [--outer 5] [--inner 5] [--repeats 1] [--seed 0]
Output:  nested AUC (mean, SD over outer folds), the flat score (best inner score on all data, optimistic),
         the number of configurations searched, and the configuration chosen in each outer fold.
'''
import argparse

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def make_pipeline():
    return Pipeline([('scaler', StandardScaler()),
                     ('select', SelectKBest(f_classif)),
                     ('clf', LogisticRegression(max_iter=5000))])


def splitter(groups, n_splits, seed):
    if groups is None:
        return StratifiedKFold(n_splits, shuffle=True, random_state=seed)
    return StratifiedGroupKFold(n_splits, shuffle=True, random_state=seed)


def nested_auc(X, y, groups, ks, Cs, outer, inner, repeats, seed):
    grid = {'select__k': ks, 'clf__C': Cs}
    scores, chosen = [], []
    for repeat in range(repeats):
        outer_cv = splitter(groups, outer, seed + 1000 + repeat)
        for train, test in outer_cv.split(X, y, groups):
            search = GridSearchCV(make_pipeline(), grid, cv=splitter(groups, inner, seed + repeat),
                                  scoring='roc_auc')
            fit_args = {} if groups is None else {'groups': groups[train]}
            search.fit(X[train], y[train], **fit_args)
            if len(set(y[test])) < 2:
                continue
            scores.append(roc_auc_score(y[test], search.predict_proba(X[test])[:, 1]))
            chosen.append(search.best_params_)
    flat = GridSearchCV(make_pipeline(), grid, cv=splitter(groups, inner, seed), scoring='roc_auc')
    flat.fit(X, y, **({} if groups is None else {'groups': groups}))
    return np.array(scores), flat.best_score_, chosen, len(ks) * len(Cs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('csv')
    parser.add_argument('--label', required=True)
    parser.add_argument('--group', help='patient/donor column; rows sharing it never split across folds')
    parser.add_argument('--drop', nargs='*', default=[], help='non-feature columns to ignore, e.g. sample_id')
    parser.add_argument('--ks', default='10,50,200', help='comma-separated feature counts to search')
    parser.add_argument('--Cs', default='0.01,0.1,1', help='comma-separated C values to search')
    parser.add_argument('--outer', type=int, default=5)
    parser.add_argument('--inner', type=int, default=5)
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()

    data = pd.read_csv(args.csv)
    non_features = [args.label, *args.drop] + ([args.group] if args.group else [])
    features = data.drop(columns=non_features).select_dtypes('number')
    y = data[args.label].to_numpy()
    groups = data[args.group].to_numpy() if args.group else None
    if groups is None:
        for column in args.drop:
            if data[column].duplicated().any():
                print(f'WARNING: column {column!r} has repeated values; if rows share a patient or donor, '
                      f'rerun with --group {column} or the estimate leaks')
    ks = sorted({min(int(k), features.shape[1]) for k in args.ks.split(',')})
    Cs = [float(c) for c in args.Cs.split(',')]
    scores, flat, chosen, n_config = nested_auc(features.to_numpy(float), y, groups, ks, Cs,
                                                args.outer, args.inner, args.repeats, args.seed)
    print(f'{args.csv}: {len(y)} rows, {features.shape[1]} features'
          + (f', {len(set(groups))} groups' if groups is not None else ''))
    print(f'Nested ROC AUC = {scores.mean():.2f} (SD {scores.std(ddof=1) if len(scores) > 1 else 0:.2f} over '
          f'{len(scores)} outer folds; {n_config} configurations searched)')
    print(f'Flat selected score = {flat:.2f}  (best inner CV score on all data; optimistic, do not report as performance)')
    for params in chosen:
        print('chosen:', params)
