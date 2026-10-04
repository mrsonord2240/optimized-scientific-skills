'''Boruta all-relevant selection keeps the whole correlated module (unlike LASSO).

Runs end-to-end on synthetic data. Five co-expressed genes (g0-g4) all carry the
signal. An all-relevant selector (Boruta) is compared with a cross-validated L1
(minimal-optimal) fit on the same data; the script prints how many module members
each keeps, because L1 may drop redundant members and absence from its list then
means "redundant", not "irrelevant".
'''
# Reference: numpy 2.5, pandas 3.0, scikit-learn 1.9, boruta 0.4.3 | Verify API if version differs

import numpy as np
import pandas as pd
from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

rng = np.random.default_rng(0)
n, p = 200, 60
latent = rng.normal(size=n)                            # shared biological signal
X = rng.normal(size=(n, p))
X[:, :5] = latent[:, None] + rng.normal(scale=0.3, size=(n, 5))   # g0-g4: one co-expressed module
y = (latent + rng.normal(scale=0.3, size=n) > 0).astype(int)
X = pd.DataFrame(X, columns=[f'g{i}' for i in range(p)])

# Boruta needs a tree estimator and numpy arrays. perc=100 uses the max shadow importance.
rf = RandomForestClassifier(n_estimators=200, n_jobs=-1, max_depth=5, random_state=42)
boruta = BorutaPy(rf, n_estimators='auto', perc=100, two_step=True, max_iter=100, random_state=42)
boruta.fit(X.values, y)                                # numpy arrays, not pandas

confirmed = list(X.columns[boruta.support_])
tentative = list(X.columns[boruta.support_weak_])
print(f'Confirmed all-relevant ({len(confirmed)}): {confirmed}')
print(f'Tentative ({len(tentative)}): {tentative}')
module = ['g0', 'g1', 'g2', 'g3', 'g4']
module_recovered = sum(g in confirmed for g in module)
print(f'\nBoruta (all-relevant) module members g0-g4 confirmed: {module_recovered}/5')

# Minimal-optimal comparison: L1 logistic, strength chosen by 5-fold CV on AUC, scaled inside the pipeline.
l1 = make_pipeline(StandardScaler(),
                   LogisticRegressionCV(l1_ratios=[1], solver='liblinear', Cs=10, cv=5,
                                        scoring='roc_auc', use_legacy_attributes=False)).fit(X, y)
kept = list(X.columns[l1[-1].coef_[0] != 0])
print(f'L1 (minimal-optimal) kept {len(kept)} genes, module members g0-g4: '
      f'{sum(g in kept for g in module)}/5')
