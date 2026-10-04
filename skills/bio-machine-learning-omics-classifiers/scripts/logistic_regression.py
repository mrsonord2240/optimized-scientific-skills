'''Elastic-net logistic fit, plus a batch-shortcut sanity check.

Runs end-to-end on synthetic data. Part 1 fits an elastic-net logistic model and reports
how many coefficients are non-zero. Part 2 demonstrates the most important omics-classifier
failure mode: when batch is confounded with the outcome, a random-split CV reaches a high AUC
by learning the batch. The batch-prediction check flags it and leave-one-batch-out
(batch_checks.py) returns the honest per-batch AUC.
'''
# Reference: numpy 2.5, pandas 3.0, scikit-learn 1.9 | Verify API if version differs

import numpy as np
from batch_checks import batch_predictability, leave_one_batch_out
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression

rng = np.random.default_rng(0)

# --- Part 1: elastic-net signature ---
X, y = make_classification(n_samples=300, n_features=500, n_informative=15, random_state=0)
X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.3, stratify=y, random_state=0)
# saga supports mixed L1/L2 (elastic net); standardize because the penalty is scale-sensitive.
pipe = Pipeline([('s', StandardScaler()),
                 ('c', LogisticRegression(solver='saga', l1_ratio=0.5,
                                          C=0.1, max_iter=5000))])
pipe.fit(X_tr, y_tr)
n_sel = int((pipe.named_steps['c'].coef_[0] != 0).sum())
print(f'Elastic-net selected {n_sel} of {X.shape[1]} features')

# --- Part 2: batch confounding has no computational rescue ---
# The features contain ZERO biological signal. Batch is confounded with the label
# (batches 0-2 mostly controls, 3-5 mostly cases) and adds a strong technical
# offset. A high CV AUC here is pure artifact: there is nothing real to learn.
n_per, n_batch = 40, 6
batch = np.repeat(np.arange(n_batch), n_per)
case_rate = np.where(batch < 3, 0.15, 0.85)            # label correlated with batch
label = (rng.random(n_per * n_batch) < case_rate).astype(int)
tech = rng.normal(batch[:, None] * 2.0, 1.0, (len(label), 100))         # batch offset, NO biology
Xb = tech

clf = Pipeline([('s', StandardScaler()), ('c', LogisticRegression(max_iter=5000))])
naive = cross_val_score(clf, Xb, label, cv=5, scoring='roc_auc')
batch_pred = batch_predictability(clf, Xb, batch)                           # can we predict the BATCH?
pooled, per_batch = leave_one_batch_out(clf, Xb, label, batch)

print(f'\nRandom-split CV AUC on data with ZERO real signal: {naive.mean():.2f}  <- batch artifact')
print(f'Batch predictability AUC (6 batches, one-vs-rest): {batch_pred:.2f}  (high = the red flag)')
print(f'Leave-one-batch-out mean per-batch AUC:            {per_batch.auc.mean():.2f}  (chance = 0.50; pooled {pooled:.2f})')
print('Confounded batch has no computational rescue (Soneson 2014) -- fix at the design stage.')
