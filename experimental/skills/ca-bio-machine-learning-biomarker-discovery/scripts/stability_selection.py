'''Stability selection with a fixed per-subsample budget q, independent of feature units.

Each subsample is standardized on its own rows, then an L1 logistic path is walked from
the strongest penalty down until more than q coefficients are non-zero; the last fit with
at most q selected is that subsample's selection (Meinshausen-Buhlmann 2010 formulation).
Regularization is therefore set by the error budget q, not by a unit-dependent C.
Identical calls (same data, q, seed) return identical frequencies: both the subsampling and the
liblinear coordinate order derive from `seed`. Unit invariance is approximate: rescaling
features changes the path's floating-point arithmetic, so a marginal feature can flip.
'''
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.svm import l1_min_c


def stability_selection(X, y, q=30, n_subsample=100, seed=0, grid_step=1.25, max_steps=60):
    X = np.asarray(X, dtype=float); y = np.asarray(y)
    n, p = X.shape
    rng = np.random.default_rng(seed)
    Z = np.zeros((n_subsample, p), dtype=int)
    for b in range(n_subsample):
        idx = rng.choice(n, size=n // 2, replace=False)           # n/2 subsampling
        Xb, yb = X[idx], y[idx]
        sd = Xb.std(axis=0); sd[sd == 0] = 1.0
        Xb = (Xb - Xb.mean(axis=0)) / sd                          # scale inside the resample
        C = l1_min_c(Xb, yb, loss='log') * grid_step
        fit_seed = int(rng.integers(2**31 - 1))                   # liblinear shuffles coordinates: seed it per subsample
        for _ in range(max_steps):
            fit = LogisticRegression(l1_ratio=1, solver='liblinear', C=C, max_iter=2000,
                                     random_state=fit_seed).fit(Xb, yb)
            mask = fit.coef_[0] != 0
            if mask.sum() > q:
                break
            Z[b] = mask
            C *= grid_step
    freq = Z.mean(axis=0)
    k = Z.sum(axis=1)
    kbar = k.mean()
    # Nogueira 2018 stability index; undefined when nothing (or everything) is selected.
    if kbar == 0 or kbar == p:
        stability = float('nan')
    else:
        stability = 1 - Z.var(axis=0, ddof=1).mean() / ((kbar / p) * (1 - kbar / p))
    return freq, stability
