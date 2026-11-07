'''LIME explanation of ONE prediction, repeated across seeds, with top-k overlap.

Inputs:  CSV, one row per sample, a binary label column, numeric feature columns.
Usage:   python scripts/lime_local.py data.csv --label label --row 0 [--positive AML] [--drop id,batch]
             [--seeds 5] [--num-features 5] [--no-discretize] [--kernel-width W] [--max-features 2000]
             [--out lime_results.tsv]
What it does: fits a 300-tree random forest on a 70% training split, then explains row `--row` of the file with
         LIME under `--seeds` different random seeds. The top-k feature sets are compared: a ranking that moves
         with the seed is not a finding. LIME is local, perturbs off the data manifold for correlated genes, and
         can be fooled by a deliberately biased model; never use it for a global ranking or as evidence a model is
         unbiased.
Output:  TSV: one row per (seed, feature) with the local weight, then summary rows (distinct top-k sets, mean
         pairwise top-k overlap). Refuses a model with test AUC under 0.6.
'''
import argparse
import sys

import numpy as np
import pandas as pd
from omics_data import load_xy
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split


def main():
    from lime.lime_tabular import LimeTabularExplainer  # imported late so --help works without it

    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('csv')
    ap.add_argument('--label', required=True)
    ap.add_argument('--row', type=int, required=True, help='0-based row of the CSV to explain')
    ap.add_argument('--positive')
    ap.add_argument('--drop', default='')
    ap.add_argument('--seeds', type=int, default=5)
    ap.add_argument('--num-features', type=int, default=5)
    ap.add_argument('--num-samples', type=int, default=5000)
    ap.add_argument('--no-discretize', action='store_true')
    ap.add_argument('--kernel-width', type=float)
    ap.add_argument('--max-features', type=int, default=2000, help='keep the highest-variance features (label-free) above this count')
    ap.add_argument('--model-seed', type=int, default=0)
    ap.add_argument('--out', default='lime_results.tsv')
    a = ap.parse_args()
    X, y, _ = load_xy(a.csv, a.label, a.positive, a.drop)
    if not 0 <= a.row < len(X):
        sys.exit(f'REFUSED: --row {a.row} is outside 0..{len(X) - 1}')
    if a.seeds < 2:
        sys.exit('REFUSED: --seeds must be at least 2; a single LIME run says nothing about its stability')
    if X.shape[1] > a.max_features:
        X = X[X.var().nlargest(a.max_features).index]
    idx = np.arange(len(X))
    tr, te = train_test_split(idx, test_size=0.3, stratify=y, random_state=a.model_seed)
    model = RandomForestClassifier(n_estimators=300, n_jobs=-1, random_state=a.model_seed).fit(X.iloc[tr], y[tr])
    auc = roc_auc_score(y[te], model.predict_proba(X.iloc[te])[:, 1])
    print(f'random forest test AUC {auc:.3f}')
    if auc < 0.6:
        sys.exit('REFUSED: the model does not generalize (test AUC under 0.6); a local explanation would explain an artifact.')
    if a.row in set(tr.tolist()):
        print('NOTE: the explained row was in the training split.')
    rows, sets = [], []
    for seed in range(a.seeds):
        explainer = LimeTabularExplainer(X.iloc[tr].to_numpy(), feature_names=list(X.columns), mode='classification',
                                         discretize_continuous=not a.no_discretize, kernel_width=a.kernel_width, random_state=seed)
        exp = explainer.explain_instance(X.iloc[a.row].to_numpy(), model.predict_proba, num_features=a.num_features, num_samples=a.num_samples)
        local = exp.local_exp[1]                               # (feature index, weight) pairs for the positive class
        sets.append(frozenset(X.columns[i] for i, _ in local))
        rows += [{'seed': seed, 'feature_rule': X.columns[i], 'weight': w} for i, w in local]
    table = pd.DataFrame(rows)
    pairs = [len(sets[i] & sets[j]) for i in range(a.seeds) for j in range(i + 1, a.seeds)]
    summary = pd.DataFrame([
        {'seed': 'summary', 'feature_rule': 'distinct_top_sets', 'weight': len(set(sets))},
        {'seed': 'summary', 'feature_rule': f'mean_pairwise_overlap_of_{a.num_features}', 'weight': float(np.mean(pairs))}])
    pd.concat([table, summary], ignore_index=True).to_csv(a.out, sep='\t', index=False)
    print(f'{len(set(sets))} distinct top-{a.num_features} sets across {a.seeds} seeds for the SAME prediction; '
          f'mean pairwise overlap {np.mean(pairs):.2f}/{a.num_features}')
    print('This is a sensitivity measurement for one prediction, not a universal claim. Keep LIME local.')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
