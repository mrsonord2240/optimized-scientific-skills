'''Predictive survival modeling: penalized Cox vs RSF, evaluated beyond the C-index.

Fits an elastic-net Cox (alpha chosen by cross-validation inside the training split
only) and a random survival forest, then evaluates on the held-out split with Uno's
IPCW C, time-dependent AUC, and integrated Brier score (IBS) versus a censoring-aware
Kaplan-Meier baseline scored with the same estimator, time grid and split.

Usage: python scripts/cox_regression.py [--data synthetic|gbsg2]
  synthetic (default): seeded simulation. gbsg2: breast-cancer cohort bundled with scikit-survival (686 patients).
'''
# Reference: scikit-survival 0.28, scikit-learn 1.9, numpy 2.5 | Verify API if version differs

import argparse

import numpy as np
from sklearn.model_selection import KFold, train_test_split
from sksurv.ensemble import RandomSurvivalForest
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored, concordance_index_ipcw, cumulative_dynamic_auc, integrated_brier_score
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.util import Surv


def load_data(name):
    '''Return X (ndarray), y (structured: bool event, float time).'''
    if name == 'gbsg2':
        from sksurv.datasets import load_gbsg2
        from sksurv.preprocessing import OneHotEncoder
        X, y = load_gbsg2()
        X = OneHotEncoder().fit_transform(X).astype(float).to_numpy()
        return X, Surv.from_arrays(event=y['cens'].astype(bool), time=y['time'].astype(float))
    rng = np.random.default_rng(0)
    n, p = 600, 40
    X = rng.normal(size=(n, p))
    beta = np.zeros(p)
    beta[:5] = 0.8                                       # 5 prognostic features
    true_time = rng.exponential(np.exp(-(X @ beta)))     # higher risk -> shorter time
    censor = rng.exponential(1.5, n)
    return X, Surv.from_arrays(event=true_time <= censor, time=np.minimum(true_time, censor))


def fit_coxnet_cv(X_train, y_train, l1_ratio=0.9, n_splits=5, seed=0):
    '''Elastic-net Cox with alpha chosen by K-fold CV (Harrell's C on the validation folds) on the TRAINING data only.

    The returned model is refit on all training data at the single selected alpha, so
    predict() and predict_survival_function() use that alpha (a default path fit would
    predict at the smallest alpha on the path, the least-penalized model).
    '''
    path = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alpha_min_ratio=0.01).fit(X_train, y_train).alphas_
    cv_c = np.zeros((n_splits, len(path)))
    for k, (i, j) in enumerate(KFold(n_splits, shuffle=True, random_state=seed).split(X_train)):
        m = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alphas=path).fit(X_train[i], y_train[i])   # one path fit per fold
        cv_c[k] = [concordance_index_censored(y_train[j]['event'], y_train[j]['time'],
                                              m.predict(X_train[j], alpha=a))[0] for a in path]
    best_alpha = path[cv_c.mean(axis=0).argmax()]
    return CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alphas=[best_alpha], fit_baseline_model=True).fit(X_train, y_train)


def km_baseline_surv(y_train, times, n_rows):
    '''Censoring-aware Kaplan-Meier survival at `times` from the training data, tiled to n_rows.'''
    t_km, s_km = kaplan_meier_estimator(y_train['event'], y_train['time'])
    s_at = np.r_[1.0, s_km][np.searchsorted(t_km, times, side='right')]   # step function; S=1 before the first time
    return np.tile(s_at, (n_rows, 1))


def main(data='synthetic'):
    X, y = load_data(data)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=0, stratify=y['event'])
    print(f'data={data}  n_train={len(y_train)}  n_test={len(y_test)}  p={X.shape[1]}  '
          f'censored in train={1 - y_train["event"].mean():.0%}')

    cox = fit_coxnet_cv(X_train, y_train)
    rsf = RandomSurvivalForest(n_estimators=300, min_samples_leaf=15, random_state=0, n_jobs=-1).fit(X_train, y_train)
    print(f'elastic-net Cox: CV-selected alpha={cox.alphas_[0]:.4f}, '
          f'{int((cox.coef_ != 0).sum())} of {X.shape[1]} coefficients nonzero')

    t_horizon = np.percentile(y_train['time'][y_train['event']], 90)           # Uno C truncation tau
    times = np.percentile(y_test['time'][y_test['event']], np.linspace(10, 80, 15))   # inside follow-up

    print(f'\n{"model":18s} {"UnoC":>6s} {"meanAUC":>8s} {"IBS":>6s}')
    for name, m in [('elastic-net Cox', cox), ('random forest', rsf)]:
        risk = m.predict(X_test)                         # higher = higher risk, NOT a probability
        c = concordance_index_ipcw(y_train, y_test, risk, tau=t_horizon)[0]
        _, mean_auc = cumulative_dynamic_auc(y_train, y_test, risk, times)
        surv = np.vstack([[fn(t) for t in times] for fn in m.predict_survival_function(X_test)])
        ibs = integrated_brier_score(y_train, y_test, surv, times)
        print(f'{name:18s} {c:6.3f} {mean_auc:8.3f} {ibs:6.3f}')

    # KM-only baseline: a model whose IBS does not beat this has no predictive value.
    km_ibs = integrated_brier_score(y_train, y_test, km_baseline_surv(y_train, times, len(y_test)), times)
    print(f'\nKaplan-Meier baseline IBS (no covariates): {km_ibs:.3f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', choices=['synthetic', 'gbsg2'], default='synthetic')
    main(ap.parse_args().data)
