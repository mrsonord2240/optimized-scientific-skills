#!/usr/bin/env python3
'''Honest ROC AUC of a FIXED pipeline (no tuning) by repeated k-fold, grouped when rows share a patient.

Every data-dependent step (univariate selection, scaling, the classifier) is refit inside each training
fold, so one CV loop is enough: nested CV is only needed when something is tuned or chosen.

Inputs:  CSV with one row per sample, a binary label column, optional sample-id and group columns, and
         numeric feature columns (everything else).
Usage:   python cv_auc.py data.csv --label label [--group patient_id] [--drop sample_id] [--k 50]
                          [--C 1.0] [--folds 5] [--repeats 10] [--seed 0]
         or import cv_auc and call repeated_cv_auc(pipe, X, y, groups) with your own Pipeline.
Output:  one line: mean AUC over repeats, its SD, and the scoring unit. That mean is the estimate.
Checked: scikit-learn 1.9.1, pandas 2.3.3, numpy 2.5.3. 60 x 5000 and 120 x 2000 inputs each run in seconds.
'''
import argparse

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def default_pipeline(k=50, C=1.0):
    # f_classif on two classes ranks features exactly as a two-sample t-test does (F = t^2)
    return Pipeline([('select', SelectKBest(f_classif, k=k)),
                     ('scale', StandardScaler()),
                     ('clf', LogisticRegression(C=C, max_iter=5000))])


def _splits(y, groups, folds, seed):
    '''Row-index splits. With groups, whole groups move together and folds are stratified by group label.'''
    rows = np.arange(len(y))
    if groups is None:
        yield from StratifiedKFold(folds, shuffle=True, random_state=seed).split(rows, y)
        return
    per_group = pd.Series(y).groupby(groups).nunique()
    if (per_group > 1).any():  # label varies within a group: best-effort stratification
        yield from StratifiedGroupKFold(folds, shuffle=True, random_state=seed).split(rows, y, groups)
        return
    names = per_group.index.to_numpy()
    group_y = pd.Series(y).groupby(groups).first().loc[names].to_numpy()
    for train_g, test_g in StratifiedKFold(folds, shuffle=True, random_state=seed).split(names, group_y):
        yield rows[np.isin(groups, names[train_g])], rows[np.isin(groups, names[test_g])]


def repeated_cv_auc(pipe, X, y, groups=None, folds=5, repeats=10, seed=0):
    '''Return one AUC per repeat, from that repeat's pooled out-of-fold scores.

    With groups whose label is constant, scores are averaged per group first, so the unit is the group.'''
    X, y = np.asarray(X, dtype=float), np.asarray(y)
    groups = None if groups is None else np.asarray(groups)
    by_group = groups is not None and (pd.Series(y).groupby(groups).nunique() == 1).all()
    aucs = []
    for repeat in range(repeats):
        scores = np.full(len(y), np.nan)
        for train, test in _splits(y, groups, folds, seed + repeat):
            model = clone(pipe).fit(X[train], y[train])
            scores[test] = model.decision_function(X[test])
        if by_group:
            frame = pd.DataFrame({'g': groups, 'y': y, 's': scores}).groupby('g').agg(y=('y', 'first'), s=('s', 'mean'))
            aucs.append(roc_auc_score(frame['y'], frame['s']))
        else:
            aucs.append(roc_auc_score(y, scores))
    return np.array(aucs), ('group' if by_group else 'sample')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('csv')
    parser.add_argument('--label', required=True)
    parser.add_argument('--group', help='column naming the patient/donor; rows sharing it never split across folds')
    parser.add_argument('--drop', nargs='*', default=[], help='non-feature columns to ignore, e.g. sample_id')
    parser.add_argument('--k', type=int, default=50)
    parser.add_argument('--C', type=float, default=1.0)
    parser.add_argument('--folds', type=int, default=5)
    parser.add_argument('--repeats', type=int, default=10)
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
    aucs, unit = repeated_cv_auc(default_pipeline(min(args.k, features.shape[1]), args.C), features.to_numpy(),
                                 y, groups, args.folds, args.repeats, args.seed)
    print(f'{args.csv}: {len(y)} rows, {features.shape[1]} features'
          + (f', {len(set(groups))} groups' if groups is not None else ''))
    print(f'ROC AUC = {aucs.mean():.2f} (SD {aucs.std(ddof=1) if len(aucs) > 1 else 0:.2f} over {args.repeats} repeats '
          f'of {args.folds}-fold; scored per {unit}; selection and scaling refit in every training fold)')
