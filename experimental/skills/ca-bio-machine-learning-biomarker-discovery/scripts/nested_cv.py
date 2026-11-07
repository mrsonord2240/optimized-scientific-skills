"""Nested CV for independent-sample numeric biomarker data."""
import argparse
import json
from pathlib import Path
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input")
p.add_argument("outcome")
p.add_argument("output_json")
args = p.parse_args()
data = pd.read_csv(args.input)
if args.outcome not in data:
    p.error("outcome column missing")
y = data.pop(args.outcome)
if y.isna().any() or y.nunique() != 2 or y.value_counts().min() < 10:
    p.error("binary outcome needs at least ten observations per class")
X = data.apply(pd.to_numeric, errors="raise")
if X.isna().any().any() or X.shape[1] < 1:
    p.error("features must be complete and numeric")
ks = [k for k in (10, 20, 50) if k <= X.shape[1]] or [X.shape[1]]
pipe = Pipeline([("scale", StandardScaler()), ("select", SelectKBest(f_classif)),
                 ("clf", LogisticRegression(max_iter=5000))])
inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
outer = StratifiedKFold(n_splits=10, shuffle=True, random_state=0)
search = GridSearchCV(pipe, {"select__k": ks}, scoring="roc_auc", cv=inner)
auc = cross_val_score(search, X, y, cv=outer, scoring="roc_auc")
Path(args.output_json).write_text(json.dumps({"fold_auc": auc.tolist(), "mean_auc": float(auc.mean()), "k_grid": ks}, indent=2) + "\n", encoding="utf-8")
