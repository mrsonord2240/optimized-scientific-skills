'''Elastic-net logistic baseline for omics data, optionally a lasso signature.

Inputs:  CSV, one row per sample, a binary label column, numeric feature columns.
Usage:   python scripts/logistic_regression.py data.csv --label label [--positive AML] [--drop id,batch]
             [--signature] [--l1-ratios 0.1,0.5,0.9] [--test-size 0.3] [--seed 0]
             [--out baseline_results.tsv] [--predictions test_predictions.csv]
Output:  metric/value TSV: non-zero coefficients, test AUC, Brier, calibration slope, chosen C and l1_ratio.
         The predictions CSV (label, prob) feeds model-validation's calibration route.
Scaling is fit inside the pipeline, so it never sees the test split. One hold-out split is a fit check, not a
performance estimate: use model-validation for that.
'''
import argparse
import sys

import numpy as np
from omics_data import load_xy
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def calibration_slope(y, p):
    q = np.clip(p, 1e-6, 1 - 1e-6)
    logit = np.log(q / (1 - q))
    return float(LogisticRegression(C=1e6, max_iter=1000).fit(logit[:, None], y).coef_[0, 0])


def fit_dense(X_tr, y_tr, l1_ratios, seed):
    clf = LogisticRegressionCV(solver='saga', l1_ratios=l1_ratios, Cs=20, cv=5, max_iter=10000,
                               scoring='neg_log_loss', use_legacy_attributes=False, random_state=seed)
    return Pipeline([('scaler', StandardScaler()), ('clf', clf)]).fit(X_tr, y_tr)


def fit_signature(X_tr, y_tr, seed):
    '''Lasso at the strongest penalty within one standard error of the best CV log loss.'''
    lasso = LogisticRegressionCV(solver='saga', l1_ratios=[1.0], Cs=20, cv=5, max_iter=10000,
                                 scoring='neg_log_loss', use_legacy_attributes=False, random_state=seed)
    Pipeline([('scaler', StandardScaler()), ('clf', lasso)]).fit(X_tr, y_tr)
    folds = lasso.scores_[:, 0, :]                       # (cv folds, Cs)
    mean, se = folds.mean(0), folds.std(0, ddof=1) / np.sqrt(folds.shape[0])
    c_1se = float(lasso.Cs_[np.flatnonzero(mean >= mean.max() - se[mean.argmax()]).min()])  # Cs_ ascends: smallest C = strongest penalty
    return Pipeline([('scaler', StandardScaler()),
                     ('clf', LogisticRegression(solver='saga', l1_ratio=1.0, C=c_1se, max_iter=10000, random_state=seed))]).fit(X_tr, y_tr), c_1se


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('csv')
    ap.add_argument('--label', required=True)
    ap.add_argument('--positive')
    ap.add_argument('--drop', default='')
    ap.add_argument('--signature', action='store_true', help='lasso at the 1-SE penalty for a short feature list')
    ap.add_argument('--l1-ratios', default='0.1,0.5,0.9')
    ap.add_argument('--test-size', type=float, default=0.3)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default='baseline_results.tsv')
    ap.add_argument('--predictions', default='test_predictions.csv')
    a = ap.parse_args()
    X, y, _ = load_xy(a.csv, a.label, a.positive, a.drop)
    if not 0.05 <= a.test_size <= 0.5:
        sys.exit('REFUSED: --test-size must be between 0.05 and 0.5')
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=a.test_size, stratify=y, random_state=a.seed)
    if min(y_tr.sum(), (1 - y_tr).sum()) < 5:
        sys.exit('REFUSED: fewer than 5 samples of one class in the training split; the 5-fold penalty search cannot run')
    if a.signature:
        model, c_used = fit_signature(X_tr, y_tr, a.seed)
        l1_used = 1.0
    else:
        model = fit_dense(X_tr, y_tr, [float(v) for v in a.l1_ratios.split(',')], a.seed)
        c_used, l1_used = float(model['clf'].C_), float(model['clf'].l1_ratio_)
    coef = model['clf'].coef_[0]
    p = model.predict_proba(X_te)[:, 1]
    rows = [('model', 'lasso signature (1-SE)' if a.signature else 'elastic-net logistic'),
            ('n_train', len(y_tr)), ('n_test', len(y_te)), ('n_features', X.shape[1]),
            ('nonzero_coefficients', int((coef != 0).sum())), ('C', c_used), ('l1_ratio', l1_used),
            ('test_auc', roc_auc_score(y_te, p)), ('test_brier', brier_score_loss(y_te, p)),
            ('test_calibration_slope', calibration_slope(y_te, p)),
            ('test_prevalence', float(y_te.mean())), ('mean_predicted', float(p.mean()))]
    with open(a.out, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('metric\tvalue\n')
        for k, v in rows:
            fh.write(f'{k}\t{v}\n')
            print(f'{k}\t{v}')
    with open(a.predictions, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('label,prob\n')
        fh.writelines(f'{int(t)},{q:.6f}\n' for t, q in zip(y_te, p))
    top = X.columns[np.argsort(-np.abs(coef))[:min(10, int((coef != 0).sum()))]].tolist()
    print('strongest coefficients:', top)
    print('wrote', a.out, 'and', a.predictions)


if __name__ == '__main__':
    main()
