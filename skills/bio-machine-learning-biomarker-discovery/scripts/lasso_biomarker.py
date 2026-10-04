'''Leakage-safe selection and stability, demonstrated on synthetic p>>n data.

Three teaching points run end-to-end on synthetic data:
1. Selecting features on the FULL matrix before CV inflates AUC on pure-noise labels;
   putting scaling and selection inside the CV fold does not.
2. A signature is reported with a stability index, not on its own.
3. The stable set is reported next to its permuted-label null count, and does not
   change when feature units change.
'''
# Reference: numpy 2.5, pandas 3.0, scikit-learn 1.9 | Verify API if version differs

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold

from stability_selection import stability_selection

rng = np.random.default_rng(0)
n, p = 80, 5000                                  # p >> n: 80 samples, 5000 'genes'
X = pd.DataFrame(rng.normal(size=(n, p)), columns=[f'g{i}' for i in range(p)])
y_noise = rng.integers(0, 2, n)                  # pure noise: no gene is truly associated
cv = StratifiedKFold(10, shuffle=True, random_state=0)

# WRONG: select top-20 on all data, then CV the classifier on those 20.
top20 = SelectKBest(f_classif, k=20).fit(X, y_noise).get_support(indices=True)
leaky = cross_val_score(LogisticRegression(max_iter=5000), X.iloc[:, top20], y_noise,
                        cv=cv, scoring='roc_auc')

# RIGHT: scaling and selection inside the Pipeline, re-fit per fold.
pipe = Pipeline([('scale', StandardScaler()),
                 ('select', SelectKBest(f_classif, k=20)),
                 ('clf', LogisticRegression(max_iter=5000))])
honest = cross_val_score(pipe, X, y_noise, cv=cv, scoring='roc_auc')

print(f'AUC on PURE NOISE, selection-before-CV (WRONG): {leaky.mean():.2f}')   # chance is 0.5
print(f'AUC on PURE NOISE, selection-inside-CV (RIGHT): {honest.mean():.2f}')

# --- Stability on a separate matrix with real signal in the first 5 genes.
# --- Moderate p after a notional pre-filter so L1 can recover the signal. ---
ns, ps, q = 120, 300, 10
Xs = pd.DataFrame(rng.normal(size=(ns, ps)), columns=[f'g{i}' for i in range(ps)])
beta = np.zeros(ps); beta[:5] = 2.5
logit = Xs.values @ beta + rng.normal(scale=1.0, size=ns)
ys = (logit > np.median(logit)).astype(int)

freq, stability = stability_selection(Xs, ys, q=q)
stable = list(Xs.columns[freq > 0.6])                           # pi_thr = 0.6
print(f'\nStable features (>60%): {stable}  [planted signal: g0..g4; a missing one was not recovered]')
print(f'Nogueira stability index: {stability:.2f}')

freq_null, _ = stability_selection(Xs, rng.permutation(ys), q=q)
print(f'Stable features under PERMUTED labels: {int((freq_null > 0.6).sum())}')

units = np.exp(rng.normal(0, 2, ps))                            # arbitrary per-gene units
freq_u, _ = stability_selection(Xs.values * units + 5, ys, q=q)
print(f'Stable features after rescaling every gene: {int((freq_u > 0.6).sum())} '
      f'(same set: {set(np.flatnonzero(freq_u > 0.6)) == set(np.flatnonzero(freq > 0.6))})')
