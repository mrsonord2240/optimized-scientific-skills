'''Class weighting, calibration and recalibration on one omics table.

Inputs:  CSV, one row per sample, a binary label column, numeric feature columns.
Usage:   python scripts/calibration_check.py data.csv --label label [--positive X] [--drop id]
             [--model logistic|rf] [--cal-size 0.25] [--test-size 0.3] [--seed 0] [--out calibration_results.tsv]
Splits:  train / calibration / test, stratified. The calibrator sees only the calibration split; every number is
         read on the test split.
Output:  one row per variant (unweighted, class_weight=balanced, balanced + recalibrated): test AUC, Brier,
         mean predicted / observed prevalence, calibration slope, then the reliability table of the unweighted
         model. Recalibration is sigmoid unless the calibration split has at least 1000 samples and 100 in the
         rarer class (isotonic). A calibrator with 3 or fewer distinct outputs is flagged degenerate: never
         report its scores.
'''
import argparse
import sys
import warnings

import numpy as np
import pandas as pd
from omics_data import load_xy
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ISOTONIC_MIN_N, ISOTONIC_MIN_EVENTS = 1000, 100


def recalibrate(model, X_cal, y_cal):
    '''Fit a calibrator for a FITTED model on data disjoint from its training and test data.'''
    y_cal = np.asarray(y_cal)
    events = int(min(y_cal.sum(), (1 - y_cal).sum()))
    isotonic = len(y_cal) >= ISOTONIC_MIN_N and events >= ISOTONIC_MIN_EVENTS
    calibrated = CalibratedClassifierCV(FrozenEstimator(model), method='isotonic' if isotonic else 'sigmoid')
    calibrated.fit(X_cal, y_cal)
    distinct = len(np.unique(calibrated.predict_proba(X_cal)[:, 1].round(6)))
    if distinct <= 3:
        warnings.warn(f'Calibrator is degenerate ({distinct} distinct probabilities on {len(y_cal)} samples, '
                      f'{events} in the rarer class); do not report scores from it.')
    return calibrated, ('isotonic' if isotonic else 'sigmoid'), events, distinct


def calibration_table(y, p, bins=5):
    '''Mean predicted vs observed risk in quantile bins of the predicted risk.'''
    y, p = np.asarray(y), np.asarray(p)
    edges = np.quantile(p, np.linspace(0, 1, bins + 1))
    idx = np.clip(np.searchsorted(edges, p, side='right') - 1, 0, bins - 1)
    return pd.DataFrame({'predicted': [p[idx == i].mean() for i in range(bins)],
                         'observed': [y[idx == i].mean() for i in range(bins)]})


def calibration_slope(y, p):
    q = np.clip(p, 1e-6, 1 - 1e-6)
    return float(LogisticRegression(C=1e6, max_iter=1000).fit(np.log(q / (1 - q))[:, None], y).coef_[0, 0])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('csv')
    ap.add_argument('--label', required=True)
    ap.add_argument('--positive')
    ap.add_argument('--drop', default='')
    ap.add_argument('--model', choices=['logistic', 'rf'], default='logistic')
    ap.add_argument('--cal-size', type=float, default=0.25)
    ap.add_argument('--test-size', type=float, default=0.3)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default='calibration_results.tsv')
    a = ap.parse_args()
    X, y, _ = load_xy(a.csv, a.label, a.positive, a.drop)
    X_dev, X_te, y_dev, y_te = train_test_split(X, y, test_size=a.test_size, stratify=y, random_state=a.seed)
    X_tr, X_cal, y_tr, y_cal = train_test_split(X_dev, y_dev, test_size=a.cal_size / (1 - a.test_size), stratify=y_dev, random_state=a.seed)
    if min(y_tr.sum(), (1 - y_tr).sum()) < 5:
        sys.exit('REFUSED: fewer than 5 samples of one class in the training split; use a larger or less imbalanced table')
    if min(y_cal.sum(), (1 - y_cal).sum()) < 3:
        sys.exit('REFUSED: fewer than 3 samples of one class in the calibration split; a calibrator fitted on it would be degenerate')

    def fit(class_weight):
        if a.model == 'rf':
            return RandomForestClassifier(n_estimators=500, min_samples_leaf=3, class_weight=class_weight, n_jobs=-1, random_state=a.seed).fit(X_tr, y_tr)
        clf = LogisticRegressionCV(Cs=10, l1_ratios=(0,), cv=min(5, int(min(y_tr.sum(), (1 - y_tr).sum()))), scoring='neg_log_loss',
                                   use_legacy_attributes=False, max_iter=5000, class_weight=class_weight)
        return Pipeline([('s', StandardScaler()), ('c', clf)]).fit(X_tr, y_tr)

    plain, balanced = fit(None), fit('balanced')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        recal, method, events, distinct = recalibrate(balanced, X_cal, y_cal)
    degenerate = any('degenerate' in str(w.message) for w in caught)
    prev = float(y_te.mean())
    rows = []
    for name, m in [('unweighted', plain), ('balanced', balanced), (f'balanced + recalibrated ({method})', recal)]:
        p = m.predict_proba(X_te)[:, 1]
        rows.append((name, roc_auc_score(y_te, p), brier_score_loss(y_te, p), float(p.mean()) / prev, calibration_slope(y_te, p)))
    p_plain = plain.predict_proba(X_te)[:, 1]
    raw_brier = rows[0][2]
    with open(a.out, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('variant\ttest_auc\ttest_brier\tmean_predicted_over_observed\tcalibration_slope\n')
        for r in rows:
            fh.write(f'{r[0]}\t{r[1]:.4f}\t{r[2]:.4f}\t{r[3]:.3f}\t{r[4]:.3f}\n')
            print(r)
        fh.write(f'#\tn_train={len(y_tr)}\tn_cal={len(y_cal)}\tn_test={len(y_te)}\tcal_events_rarer_class={events}\tcal_distinct_outputs={distinct}\n')
    print(f'test prevalence {prev:.3f}; calibration split {len(y_cal)} samples, {events} in the rarer class, method {method}, {distinct} distinct outputs')
    print('reliability (unweighted model, test):\n' + calibration_table(y_te, p_plain).round(3).to_string(index=False))
    if degenerate:
        print('WARNING: the recalibrated row comes from a degenerate calibrator; do not report it.')
    if rows[2][2] > rows[1][2]:
        print('NOTE: recalibration did not improve the Brier of the balanced model; compare against the uncalibrated Brier before adopting it.')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
