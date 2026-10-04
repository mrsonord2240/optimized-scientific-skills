'''Linear vs RF vs XGBoost on omics-like data: discrimination AND calibration.

Runs end-to-end on synthetic data with a purely LINEAR p>>n signal (no interactions).
Splits: train (all models), validation (XGBoost early stopping only), test (every reported
number). Prints AUC and Brier on the test split plus each model's probability spread, so
the calibration direction is read from the numbers, not assumed.
'''
# Reference: numpy 2.5, scikit-learn 1.9, xgboost 3.4 | Verify API if version differs

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss
from xgboost import XGBClassifier

rng = np.random.default_rng(0)
# n=600: at n=300 XGBoost stops at round 0 on this seed (a null model). Of learning rates 0.03/0.1/0.3 at n=600,
# 0.03 had the lowest validation logloss (0.577 vs 0.601 and 0.632); chosen on train/validation data only.
n, p, k = 600, 1500, 25                                # p >> n; 25 truly informative, linear
X = rng.normal(size=(n, p))
beta = np.zeros(p); beta[:k] = rng.normal(scale=1.2, size=k)
logit = X @ beta
y = (rng.random(n) < 1.0 / (1.0 + np.exp(-logit))).astype(int)          # linear generative model
X_dev, X_te, y_dev, y_te = train_test_split(X, y, test_size=0.3, stratify=y, random_state=0)
X_tr, X_val, y_tr, y_val = train_test_split(X_dev, y_dev, test_size=0.3, stratify=y_dev, random_state=0)

models = {
    # Genuine elastic-net (saga + l1_ratio); default-L2 leaves linear signal on the table.
    'elastic-net logistic': Pipeline([('s', StandardScaler()),
                                       ('c', LogisticRegressionCV(solver='saga', l1_ratios=[0.5], Cs=10, cv=5, max_iter=5000,
                                                                  scoring='neg_log_loss', use_legacy_attributes=False))]),
    'random forest': RandomForestClassifier(n_estimators=400, min_samples_leaf=3, random_state=0),
    # Early stopping watches the validation split only; scores below come from the untouched test split.
    'xgboost': XGBClassifier(n_estimators=2000, learning_rate=0.03, max_depth=4, subsample=0.8, colsample_bytree=0.5,
                             early_stopping_rounds=50, eval_metric='logloss', random_state=0),
}

print(f'{"model":24s} {"AUC":>6s} {"Brier":>7s}')
for name, m in models.items():
    m.fit(X_tr, y_tr, eval_set=[(X_val, y_val)], verbose=False) if name == 'xgboost' else m.fit(X_tr, y_tr)
    pr = m.predict_proba(X_te)[:, 1]
    # AUC (ranking) can be similar across models while Brier (calibration) differs.
    print(f'{name:24s} {roc_auc_score(y_te, pr):6.3f} {brier_score_loss(y_te, pr):7.3f}')

xgb = models['xgboost']
print(f'\nxgboost best round {xgb.best_iteration} of 2000 (validation logloss {xgb.best_score:.3f} on {len(y_val)} samples)')
if xgb.best_iteration == 0:
    print('!!! WARNING: best round is 0 -- validation loss never beat the first tree, so the xgboost row is a null model.\n'
          '!!! Do not read its AUC, Brier or probability spread; use more samples, or choose rounds by cross-validation.')
print('\nProbability spread on the test split (a wide spread is not by itself miscalibration; check a reliability curve):')
for name, m in models.items():
    pr = m.predict_proba(X_te)[:, 1]
    print(f'{name:24s} range [{pr.min():.2f}, {pr.max():.2f}]  sd {pr.std():.2f}  share outside [0.1, 0.9]: {np.mean((pr < 0.1) | (pr > 0.9)):.2f}')
