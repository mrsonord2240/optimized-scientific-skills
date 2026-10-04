'''Calibration helpers and two measured demos (class weights, isotonic vs sigmoid).

recalibrate: recalibrate a fitted model on a disjoint calibration set; sigmoid unless the
  calibration set is large enough for isotonic; warns when the result is degenerate.
calibration_table: mean predicted vs observed risk per predicted-risk quantile bin.

Run `python calibration_check.py`. Synthetic data with a known logistic generative model, so
observed risk is the truth. Demo 1: class_weight='balanced' moves predicted risk away from the
observed prevalence; recalibrating on held-out data restores it. Demo 2: isotonic vs sigmoid
calibration as the calibration set grows.
'''
# Reference: numpy 2.5, pandas 3.0, scikit-learn 1.9 | Verify API if version differs

import warnings

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import brier_score_loss, roc_auc_score

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
    return calibrated


def calibration_table(y, p, bins=5):
    '''Mean predicted vs observed risk in quantile bins of the predicted risk.'''
    y, p = np.asarray(y), np.asarray(p)
    edges = np.quantile(p, np.linspace(0, 1, bins + 1))
    idx = np.clip(np.searchsorted(edges, p, side='right') - 1, 0, bins - 1)
    return pd.DataFrame({'predicted': [p[idx == i].mean() for i in range(bins)],
                         'observed': [y[idx == i].mean() for i in range(bins)]})


def draw(n, p, k, intercept, seed):
    rng = np.random.default_rng(seed)
    beta = np.zeros(p)
    beta[:k] = np.random.default_rng(123).normal(scale=0.6, size=k)
    X = rng.normal(size=(n, p))
    return X, (rng.random(n) < 1 / (1 + np.exp(-(X @ beta + intercept)))).astype(int)


def calibration_slope(y, p):
    logit = np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))
    return LogisticRegression(C=1e6, max_iter=1000).fit(logit[:, None], y).coef_[0, 0]


def demo_class_weights():
    print('Demo 1: 8% prevalence, L2 logistic (C by CV on log loss), 1500 train / 1000 calibration / 20000 test, mean of 5 seeds')
    p, k, intercept = 50, 10, -3.0
    rows = []
    for seed in range(5):
        Xtr, ytr = draw(1500, p, k, intercept, seed * 10 + 1)
        Xcal, ycal = draw(1000, p, k, intercept, seed * 10 + 2)
        Xte, yte = draw(20000, p, k, intercept, seed * 10 + 3)
        def fit(class_weight):
            return LogisticRegressionCV(Cs=10, l1_ratios=(0,), cv=5, scoring='neg_log_loss', use_legacy_attributes=False,
                                        max_iter=5000, class_weight=class_weight).fit(Xtr, ytr)
        plain, balanced = fit(None), fit('balanced')
        recal = recalibrate(balanced, Xcal, ycal)
        for name, m in [('unweighted', plain), ('balanced', balanced), ('balanced + recalibrated', recal)]:
            pr = m.predict_proba(Xte)[:, 1]
            rows.append((name, roc_auc_score(yte, pr), brier_score_loss(yte, pr), pr.mean() / yte.mean(), calibration_slope(yte, pr)))
            if seed == 0 and name == 'unweighted':
                print('Reliability, unweighted (seed 0):\n' + calibration_table(yte, pr).round(3).to_string(index=False))
    df = pd.DataFrame(rows, columns=['model', 'AUC', 'Brier', 'mean predicted / observed', 'calibration slope'])
    print(df.groupby('model', sort=False).mean().round(3).to_string())


def demo_isotonic_vs_sigmoid():
    print('\nDemo 2: Brier of RF recalibrated on n_cal samples (30% prevalence, 10000 test, mean of 20 repeats)')
    p, k, intercept = 50, 10, -1.0
    rows = []
    warnings.filterwarnings('ignore', category=UserWarning)     # sklearn notes small classes per CV fold at n=20
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    for rep in range(20):
        Xtr, ytr = draw(300, p, k, intercept, rep * 10 + 1)
        Xpool, ypool = draw(1000, p, k, intercept, rep * 10 + 2)
        Xte, yte = draw(10000, p, k, intercept, rep * 10 + 3)
        rf = RandomForestClassifier(n_estimators=200, min_samples_leaf=3, random_state=0, n_jobs=-1).fit(Xtr, ytr)
        for n in (20, 100, 1000):
            row = {'n_cal': n, 'raw': brier_score_loss(yte, rf.predict_proba(Xte)[:, 1])}
            for method in ('sigmoid', 'isotonic'):
                cal = CalibratedClassifierCV(FrozenEstimator(rf), method=method).fit(Xpool[:n], ypool[:n])
                row[method] = brier_score_loss(yte, cal.predict_proba(Xte)[:, 1])
            rows.append(row)
    print(pd.DataFrame(rows).groupby('n_cal').mean().round(4).to_string())

    print('\nDegenerate case: a separable problem with 20 calibration samples')
    Xtr, ytr = draw(300, 5, 5, 0.0, 7)
    Xtr = Xtr + ytr[:, None] * 3.0
    Xc, yc = draw(20, 5, 5, 0.0, 8)
    Xc = Xc + yc[:, None] * 3.0
    rf = RandomForestClassifier(n_estimators=100, random_state=0).fit(Xtr, ytr)
    iso = CalibratedClassifierCV(FrozenEstimator(rf), method='isotonic').fit(Xc, yc)
    pr = iso.predict_proba(Xc)[:, 1]
    print(f'isotonic: {len(np.unique(pr.round(6)))} distinct probabilities; Brier on its own calibration data '
          f'{brier_score_loss(yc, pr):.3f} (not a result)')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        recalibrate(rf, Xc, yc)
    print(f'recalibrate(): chose sigmoid at n=20; warnings raised: {len(caught)}')


if __name__ == '__main__':
    demo_class_weights()
    demo_isotonic_vs_sigmoid()
