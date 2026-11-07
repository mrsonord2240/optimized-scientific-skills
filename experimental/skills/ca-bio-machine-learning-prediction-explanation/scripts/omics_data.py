'''Shared loader for the explanation scripts: one CSV, one label column, numeric feature columns.

load_xy refuses what would silently corrupt a fit: a non-binary label, non-numeric or missing feature values,
duplicate sample rows are not checked. Columns named in `drop` (ids, batch labels) are removed from the features.
'''
import sys

import numpy as np
import pandas as pd


def load_xy(path, label, positive=None, drop='', batch_col=None, allow_nan=False):
    df = pd.read_csv(path)
    if label not in df.columns:
        sys.exit(f'REFUSED: label column {label!r} not in {path}; columns start {list(df.columns)[:5]}')
    drop_cols = [c for c in drop.split(',') if c]
    absent = [c for c in drop_cols + ([batch_col] if batch_col else []) if c not in df.columns]
    if absent:
        sys.exit(f'REFUSED: columns not in {path}: {absent}')
    classes = sorted(df[label].dropna().unique().tolist(), key=str)
    if len(classes) != 2:
        sys.exit(f'REFUSED: {label!r} must have exactly 2 classes, found {len(classes)}: {classes[:6]}')
    if positive is None:
        if set(classes) != {0, 1}:
            sys.exit(f'REFUSED: label values {classes} are not 0/1; name the positive class with --positive')
        positive = 1
    elif str(positive) not in [str(c) for c in classes]:
        sys.exit(f'REFUSED: --positive {positive!r} is not one of {classes}')
    y = np.array([1 if str(v) == str(positive) else 0 for v in df[label]])
    batch = df[batch_col].to_numpy() if batch_col else None
    X = df.drop(columns=[label] + drop_cols + ([batch_col] if batch_col else []))
    bad = [c for c in X.columns if not pd.api.types.is_numeric_dtype(X[c])]
    if bad:
        sys.exit(f'REFUSED: non-numeric feature columns {bad[:5]}; name ids or labels with --drop')
    if not allow_nan and X.isna().any().any():
        sys.exit('REFUSED: missing values in the features; this script does not impute (impute inside the training fold, or use XGBoost)')
    if min(y.sum(), (1 - y).sum()) < 5:
        sys.exit(f'REFUSED: the rarer class has {min(y.sum(), (1 - y).sum())} samples; fewer than 5 cannot be split and scored')
    print(f'{path}: {len(y)} samples, {X.shape[1]} features, positive class {positive!r} = {int(y.sum())} samples')
    return X, y, batch
