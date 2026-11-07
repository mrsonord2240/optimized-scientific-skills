'''Elastic-net logistic vs random forest vs XGBoost on one omics table: discrimination and calibration.

Inputs:  CSV, one row per sample, a binary label column, numeric feature columns (NaN allowed for XGBoost only).
Usage:   python scripts/rf_xgboost_classifier.py data.csv --label label [--positive AML] [--drop id]
             [--xgb-rounds early|cv] [--seed 0] [--out trees_results.tsv]
Splits:  train (all models), validation (XGBoost early stopping only), test (every reported number).
         --xgb-rounds cv picks the round count by 5-fold CV on train+validation instead of early stopping;
         use it when the validation split is small (under about 50 samples).
Output:  one row per model: test AUC, Brier, probability range and spread, best round (XGBoost).
         A best round of 0 marks XGBoost a null model: its row must not be read.
'''
import argparse
import sys

import numpy as np
from omics_data import load_xy
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('csv')
    ap.add_argument('--label', required=True)
    ap.add_argument('--positive')
    ap.add_argument('--drop', default='')
    ap.add_argument('--xgb-rounds', choices=['early', 'cv'], default='early')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default='trees_results.tsv')
    a = ap.parse_args()
    try:
        from xgboost import XGBClassifier
    except ImportError:
        sys.exit('REFUSED: xgboost is not installed in this environment; run in one that has it')
    X, y, _ = load_xy(a.csv, a.label, a.positive, a.drop, allow_nan=True)
    X_dev, X_te, y_dev, y_te = train_test_split(X, y, test_size=0.3, stratify=y, random_state=a.seed)
    X_tr, X_val, y_tr, y_val = train_test_split(X_dev, y_dev, test_size=0.3, stratify=y_dev, random_state=a.seed)
    if min(y_tr.sum(), (1 - y_tr).sum()) < 5:
        sys.exit('REFUSED: fewer than 5 samples of one class in the training split; use a larger or less imbalanced table')
    has_nan = bool(X.isna().any().any())
    # scale_pos_weight and class_weight are left out on purpose: they shift predicted risk (see the calibration route).
    models = {
        'rf': RandomForestClassifier(n_estimators=500, max_features='sqrt', min_samples_leaf=3, n_jobs=-1, random_state=a.seed),
    }
    if not has_nan:
        models = {'elastic-net logistic': Pipeline([('s', StandardScaler()), ('c', LogisticRegressionCV(
            solver='saga', l1_ratios=[0.5], Cs=10, cv=5, max_iter=5000, scoring='neg_log_loss',
            use_legacy_attributes=False, random_state=a.seed))]), **models}
    else:
        print('NaN in the features: the logistic and random forest rows are skipped (they need complete data); XGBoost handles NaN natively.')
        models = {}
    rows, best = [], None
    for name, m in models.items():
        m.fit(X_tr, y_tr)
        rows.append((name, m.predict_proba(X_te)[:, 1], ''))
    if a.xgb_rounds == 'early':
        xgb = XGBClassifier(n_estimators=2000, learning_rate=0.03, max_depth=4, subsample=0.8, colsample_bytree=0.5,
                            reg_lambda=1.0, early_stopping_rounds=50, eval_metric='logloss', n_jobs=-1, random_state=a.seed)
        xgb.fit(X_tr, y_tr, eval_set=[(X_val, y_val)], verbose=False)
        best = int(xgb.best_iteration)
        note = f'best_round={best}; validation n={len(y_val)}'
        if best == 0:
            note += '; NULL MODEL: validation loss never beat the first tree, do not read this row'
        if len(y_val) < 50:
            note += '; small validation set stops early and flatters itself, consider --xgb-rounds cv'
        rows.append(('xgboost (early stopping)', xgb.predict_proba(X_te)[:, 1], note))
    else:
        search = GridSearchCV(XGBClassifier(learning_rate=0.03, max_depth=4, subsample=0.8, colsample_bytree=0.5,
                                            n_jobs=-1, random_state=a.seed),
                              {'n_estimators': [100, 200, 400]}, cv=5, scoring='roc_auc').fit(X_dev, y_dev)
        rows.append(('xgboost (rounds by CV)', search.predict_proba(X_te)[:, 1], f'n_estimators={search.best_params_["n_estimators"]}'))
    with open(a.out, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('model\ttest_auc\ttest_brier\tprob_min\tprob_max\tprob_sd\tshare_outside_0.1_0.9\tnote\n')
        for name, p, note in rows:
            line = (name, roc_auc_score(y_te, p), brier_score_loss(y_te, p), p.min(), p.max(), p.std(),
                    float(np.mean((p < 0.1) | (p > 0.9))), note)
            fh.write('\t'.join(f'{v:.4f}' if isinstance(v, (float, np.floating)) else str(v) for v in line) + '\n')
            print(line)
    print(f'test n={len(y_te)}; a wide probability spread is not by itself miscalibration: check a reliability curve.')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
