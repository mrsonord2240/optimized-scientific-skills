'''TreeSHAP attributions for a random forest on an omics table, aggregated over co-expression modules.

Inputs:  CSV, one row per sample, a binary label column, numeric feature columns.
Usage:   python scripts/shap_tree.py data.csv --label label [--positive AML] [--drop id,batch]
             [--mode interventional|tree_path_dependent|both] [--background-n 100] [--n-explain 200]
             [--max-features 2000] [--cut 0.7] [--seed 0] [--out shap_results.tsv]
What it does: splits train/test (70/30, stratified), fits a 300-tree random forest on train, refuses to explain a
         model with test AUC under 0.6, then runs TreeExplainer with `feature_perturbation` set EXPLICITLY (never
         'auto', which flips the estimand on whether a background is given). Interventional uses a background
         sampled from the training rows; the values are in probability units for a forest.
         Modules: |Spearman| correlation on the training rows, average linkage, cut at distance (1 - rho) = --cut.
Output:  TSV, one row per feature (level=gene) and per module (level=module): mean |SHAP| per mode, module id,
         members. Rank the modules, not genes inside a module: within-module order is not a finding.
Attribution describes the model, not biology, and a ranking is not a validated biomarker panel.
'''
import argparse
import sys

import numpy as np
import pandas as pd
from omics_data import load_xy
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split


def positive_class_values(vals):
    '''shap returns a list per class, or an (n, p, classes) array, or (n, p), depending on version.'''
    if isinstance(vals, list):
        return np.asarray(vals[1])
    vals = np.asarray(vals)
    return vals[..., 1] if vals.ndim == 3 else vals


def module_map(X_train, cut):
    abs_rho = X_train.corr(method='spearman').abs().fillna(0.0).to_numpy(copy=True)  # pandas 3 may expose a read-only view
    np.fill_diagonal(abs_rho, 1.0)
    tree = linkage(squareform(1.0 - abs_rho, checks=False), method='average')
    return dict(zip(X_train.columns, fcluster(tree, t=cut, criterion='distance')))


def main():
    import shap  # imported late so --help works without it

    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('csv')
    ap.add_argument('--label', required=True)
    ap.add_argument('--positive')
    ap.add_argument('--drop', default='')
    ap.add_argument('--mode', choices=['interventional', 'tree_path_dependent', 'both'], default='interventional')
    ap.add_argument('--background-n', type=int, default=100)
    ap.add_argument('--n-explain', type=int, default=200)
    ap.add_argument('--max-features', type=int, default=2000, help='keep the highest-variance features (unsupervised) above this count')
    ap.add_argument('--cut', type=float, default=0.7, help='module cut as 1 - |rho|; an analysis choice, report it')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default='shap_results.tsv')
    a = ap.parse_args()
    X, y, _ = load_xy(a.csv, a.label, a.positive, a.drop)
    if X.shape[1] > a.max_features:
        keep = X.var().nlargest(a.max_features).index
        print(f'{X.shape[1]} features: keeping the {a.max_features} highest-variance ones (label-free) for the module map and explainer')
        X = X[keep]
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.3, stratify=y, random_state=a.seed)
    model = RandomForestClassifier(n_estimators=300, min_samples_leaf=1, n_jobs=-1, random_state=a.seed).fit(X_tr, y_tr)
    auc = roc_auc_score(y_te, model.predict_proba(X_te)[:, 1])
    print(f'random forest test AUC {auc:.3f} (n_test={len(y_te)})')
    if auc < 0.6:
        sys.exit('REFUSED: the model does not generalize (test AUC under 0.6); attributions would explain an artifact. Fix the model first.')
    explain = X_te.iloc[:a.n_explain]
    modes = ['interventional', 'tree_path_dependent'] if a.mode == 'both' else [a.mode]
    gene_tab = pd.DataFrame(index=X.columns)
    for mode in modes:
        if mode == 'interventional':
            if len(X_tr) < 5:
                sys.exit('REFUSED: too few training rows for a background')
            background = shap.utils.sample(X_tr, min(a.background_n, len(X_tr)), random_state=a.seed)
            expl = shap.TreeExplainer(model, data=background, feature_perturbation='interventional')
        else:
            expl = shap.TreeExplainer(model, feature_perturbation='tree_path_dependent')
        try:
            vals = positive_class_values(expl.shap_values(explain))
        except Exception as exc:  # additivity failure: fix the configuration, do not disable the check
            sys.exit(f'REFUSED: SHAP failed in mode {mode}: {exc}. If this is a local-accuracy (additivity) error, check raw margin versus probability; do not disable the check.')
        gene_tab[f'mean_abs_shap_{mode}'] = np.abs(vals).mean(axis=0)
    modules = module_map(X_tr, a.cut)
    gene_tab['module'] = pd.Series(modules)
    genes = gene_tab.reset_index(names='id').assign(level='gene')
    value_cols = [c for c in gene_tab.columns if c.startswith('mean_abs_shap_')]
    mod = gene_tab.groupby('module')[value_cols].sum()
    mod['members'] = gene_tab.groupby('module').size()
    mod = mod.reset_index(names='id').assign(level='module')
    out = pd.concat([mod.sort_values(value_cols[0], ascending=False), genes.sort_values(value_cols[0], ascending=False)], ignore_index=True)
    out = out[['level', 'id', 'module', 'members'] + value_cols]
    out.round(6).to_csv(a.out, sep='\t', index=False)
    print(f'{len(set(modules.values()))} modules at cut {a.cut} (1 - |Spearman|, average linkage); top modules:')
    print(mod.sort_values(value_cols[0], ascending=False).head(8).round(4).to_string(index=False))
    if len(modes) == 2:
        k = 20
        top = [set(gene_tab[c].nlargest(k).index) for c in value_cols]
        print(f'top-{k} gene overlap between the two modes: {len(top[0] & top[1])}/{k} (the module totals are the stable quantity)')
    print('scale: probability of the positive class; attribution explains the model, not biology.')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
