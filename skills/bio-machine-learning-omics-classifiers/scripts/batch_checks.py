'''Batch-shortcut checks that work for any number of batches (>= 2).

batch_predictability: can the features predict the batch? (high AUC = shortcut risk)
leave_one_batch_out: train on all batches but one, predict the held-out batch.
  The per-batch AUCs (mean of the defined ones) are the honest estimate. The pooled AUC also
  ranks samples across batches, so it still rewards a batch-level score shift that tracks the
  batch's case rate; use it only when batches are too one-sided for per-batch AUCs.

Run `python batch_checks.py` for a self-test on synthetic data with 2, 3 and 6 batches.
'''
# Reference: numpy 2.5, pandas 3.0, scikit-learn 1.9 | Verify API if version differs

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut, StratifiedKFold, cross_val_score


def batch_predictability(model, X, batch):
    '''Mean cross-validated one-vs-rest AUC for predicting the batch (works for 2 or more batches).'''
    n_splits = min(5, int(pd.Series(batch).value_counts().min()))
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)
    return cross_val_score(model, X, batch, cv=cv, scoring='roc_auc_ovr').mean()


def leave_one_batch_out(model, X, y, batch):
    '''Return (pooled AUC, per-batch table). Each batch is predicted by a model that never saw it.

    A held-out batch with one class has no AUC; it is reported as NaN in the table but its
    predictions still enter the pooled AUC. A fold whose training batches hold one class is
    skipped and listed in the table, never silently averaged in.
    '''
    X, y, batch = np.asarray(X), np.asarray(y), np.asarray(batch)
    prob = np.full(len(y), np.nan)
    rows = []
    for train, test in LeaveOneGroupOut().split(X, y, batch):
        held = batch[test][0]
        if len(np.unique(y[train])) < 2:
            rows.append((held, len(test), np.nan, 'skipped: training batches hold one class'))
            continue
        prob[test] = clone(model).fit(X[train], y[train]).predict_proba(X[test])[:, 1]
        one_class = len(np.unique(y[test])) < 2
        auc = np.nan if one_class else roc_auc_score(y[test], prob[test])
        rows.append((held, len(test), auc, 'AUC undefined: one class in this batch' if one_class else ''))
    ok = ~np.isnan(prob)
    pooled = roc_auc_score(y[ok], prob[ok]) if ok.any() and len(np.unique(y[ok])) == 2 else np.nan
    return pooled, pd.DataFrame(rows, columns=['held_out_batch', 'n', 'auc', 'note'])


if __name__ == '__main__':
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import cross_val_score as cvs

    pipe = Pipeline([('s', StandardScaler()), ('c', LogisticRegression(max_iter=5000))])
    rng = np.random.default_rng(0)
    for n_batch in (2, 3, 6):
        # No biology: features carry only a batch offset; label is confounded with batch.
        batch = np.repeat(np.arange(n_batch), 40)
        label = (rng.random(len(batch)) < np.where(batch < n_batch / 2, 0.15, 0.85)).astype(int)
        X = rng.normal(batch[:, None] * 2.0, 1.0, (len(batch), 100))
        random_auc = cvs(pipe, X, label, cv=5, scoring='roc_auc').mean()
        pooled, table = leave_one_batch_out(pipe, X, label, batch)
        print(f'{n_batch} batches: batch predictability {batch_predictability(pipe, X, batch):.2f}, '
              f'random-split label AUC {random_auc:.2f}, leave-one-batch-out mean per-batch AUC '
              f'{table.auc.mean():.2f} (pooled {pooled:.2f})')
    # A held-out batch that contains a single class.
    batch = np.repeat([0, 1, 2], 30)
    label = np.r_[np.zeros(30), np.tile([0, 1], 15), np.ones(30)].astype(int)
    X = rng.normal(size=(90, 20)) + label[:, None] * 0.8
    pooled, table = leave_one_batch_out(pipe, X, label, batch)
    print(f'\none-class batches (0 = all controls, 2 = all cases): pooled AUC {pooled:.2f}')
    print(table.to_string(index=False))
