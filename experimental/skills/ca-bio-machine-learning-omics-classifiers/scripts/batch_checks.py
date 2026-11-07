'''Batch-shortcut checks for any number of batches (2 or more).

Inputs:  CSV with a binary label column, a batch column and numeric feature columns.
Usage:   python scripts/batch_checks.py data.csv --label label --batch-col batch [--positive X] [--drop id]
             [--out batch_results.tsv]
Output:  one TSV, `table` column names each block: summary (batch predictability AUC, label-vs-batch chi-square p,
         random-split label AUC, leave-one-batch-out mean per-batch AUC and pooled AUC) and per_batch (held-out
         batch, n, auc, note).
Read the per-batch mean, not the pooled AUC: pooling also ranks samples across batches, so a batch-level score
shift that tracks each batch's case rate still scores. A held-out batch with one class has no AUC (NaN, noted); a
fold whose training batches hold one class is skipped and listed.
'''
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


def main():
    import argparse
    import sys

    from omics_data import load_xy
    from scipy.stats import chi2_contingency
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    ap = argparse.ArgumentParser(description='Batch-shortcut checks')
    ap.add_argument('csv')
    ap.add_argument('--label', required=True)
    ap.add_argument('--batch-col', required=True)
    ap.add_argument('--positive')
    ap.add_argument('--drop', default='')
    ap.add_argument('--out', default='batch_results.tsv')
    a = ap.parse_args()
    X, y, batch = load_xy(a.csv, a.label, a.positive, a.drop, batch_col=a.batch_col)
    sizes = pd.Series(batch).value_counts()
    if len(sizes) < 2:
        sys.exit('REFUSED: one batch only; there is nothing to leave out or predict')
    if sizes.min() < 2:
        sys.exit(f'REFUSED: batches with fewer than 2 samples: {sizes[sizes < 2].index.tolist()}')
    probe = Pipeline([('scaler', StandardScaler()), ('clf', LogisticRegression(max_iter=5000))])
    pred_auc = batch_predictability(probe, X, batch)
    ct = pd.crosstab(y, batch)
    p_assoc = float(chi2_contingency(ct)[1])
    rand_auc = float(cross_val_score(probe, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=0), scoring='roc_auc').mean())
    pooled, per_batch = leave_one_batch_out(probe, X, y, batch)
    summary = pd.DataFrame([
        ('batch_predictability_auc', pred_auc), ('label_vs_batch_chi2_p', p_assoc), ('random_split_label_auc', rand_auc),
        ('lobo_mean_per_batch_auc', float(per_batch.auc.mean())), ('lobo_pooled_auc', pooled),
        ('n_batches', len(sizes))], columns=['metric', 'value'])
    out = pd.concat([summary.assign(table='summary'), per_batch.assign(table='per_batch')], ignore_index=True)
    out = out[['table', 'metric', 'value', 'held_out_batch', 'n', 'auc', 'note']]
    out.round(4).to_csv(a.out, sep='\t', index=False)
    print(summary.round(4).to_string(index=False))
    print(per_batch.round(3).to_string(index=False))
    if pred_auc > 0.9 and p_assoc < 0.05:
        print('RED FLAG: batch is predictable from the features and associated with the label. A high random-split AUC '
              'may be a batch artifact; judge by the leave-one-batch-out per-batch mean. No correction rescues a confounded design.')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
