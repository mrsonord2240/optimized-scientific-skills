"""Reference solvers for i1_ml_leakage.

`correct`  : repeated cross-validation with feature selection refit inside every training
             fold (A: stratified folds; B: folds that keep each patient's samples together).
`trap_*`   : shortcuts. Each must FAIL the frozen graders.
`variant_*`: legitimate alternative estimators. Each must PASS (tolerance sanity check).

Usage:  python solvers.py        (needs data/ and eval.json from generate.py)
"""
import sys
import warnings
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import (LeaveOneGroupOut, LeaveOneOut, StratifiedGroupKFold,
                                     StratifiedKFold)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from grading import grade  # noqa: E402

warnings.filterwarnings("ignore")
REPEATS = 10

_A = pd.read_csv(HERE / "data" / "dataset_a.csv")
_B = pd.read_csv(HERE / "data" / "dataset_b.csv")
XA, yA = _A.iloc[:, 2:].to_numpy(), _A["label"].to_numpy()
XB, yB, gB = _B.iloc[:, 3:].to_numpy(), _B["label"].to_numpy(), _B["patient_id"].to_numpy()


def pipe(C=1.0, scale=True):
    steps = [SelectKBest(f_classif, k=50)] + ([StandardScaler()] if scale else []) + \
            [LogisticRegression(C=C, max_iter=5000)]
    return make_pipeline(*steps)


def _pooled_auc(X, y, splits, C=1.0, scale=True, preselect=False):
    """Out-of-fold pooled AUC. preselect=True selects the 50 features on ALL rows first."""
    if preselect:
        keep = SelectKBest(f_classif, k=50).fit(X, y).get_support()
        X = X[:, keep]
        mk = lambda: make_pipeline(*([StandardScaler()] if scale else []), LogisticRegression(C=C, max_iter=5000))
    else:
        mk = lambda: pipe(C, scale)
    pred = np.zeros(len(y))
    for tr, te in splits:
        pred[te] = mk().fit(X[tr], y[tr]).decision_function(X[te])
    return roc_auc_score(y, pred)


def _rep(X, y, groups, k, repeats=REPEATS, **kw):
    out = []
    for r in range(repeats):
        if groups is None:
            sp = StratifiedKFold(k, shuffle=True, random_state=r).split(X, y)
        else:
            sp = StratifiedGroupKFold(k, shuffle=True, random_state=r).split(X, y, groups)
        out.append(_pooled_auc(X, y, sp, **kw))
    return float(np.mean(out))


@lru_cache(None)
def comp():
    """All component estimates, computed once."""
    return dict(
        a_honest=_rep(XA, yA, None, 5),
        a_leaky=_rep(XA, yA, None, 5, preselect=True),
        b_honest=_rep(XB, yB, gB, 5),
        b_random=_rep(XB, yB, None, 5),                         # selection inside folds, patients split
        b_sel_outside=_rep(XB, yB, gB, 5, preselect=True),      # patients kept together, selection outside
        b_both=_rep(XB, yB, None, 5, preselect=True),
    )


def ans(a, b):
    return {"auc_dataset_a": round(float(a), 2), "auc_dataset_b": round(float(b), 2)}


# ---- reference and traps -----------------------------------------------------------
def correct():
    c = comp()
    return ans(c["a_honest"], c["b_honest"])


def trap_select_first_A():
    c = comp()
    return ans(c["a_leaky"], c["b_honest"])


def trap_random_split_B():
    c = comp()
    return ans(c["a_honest"], c["b_random"])


def trap_both():
    c = comp()
    return ans(c["a_leaky"], c["b_both"])


def trap_chance():
    return ans(0.5, 0.5)


def trap_select_first_B_grouped():
    c = comp()
    return ans(c["a_honest"], c["b_sel_outside"])


def trap_resubstitution():
    """Fit on everything, score on the same rows."""
    a = pipe().fit(XA, yA)
    b = pipe().fit(XB, yB)
    return ans(roc_auc_score(yA, a.decision_function(XA)), roc_auc_score(yB, b.decision_function(XB)))


# ---- legitimate variants (should pass) -------------------------------------------
def variant_10fold():
    return ans(_rep(XA, yA, None, 10, repeats=3), _rep(XB, yB, gB, 10, repeats=3))


def variant_3fold():
    return ans(_rep(XA, yA, None, 3, repeats=5), _rep(XB, yB, gB, 3, repeats=5))


def variant_leave_one_out():
    return ans(_pooled_auc(XA, yA, LeaveOneOut().split(XA, yA)),
               _pooled_auc(XB, yB, LeaveOneGroupOut().split(XB, yB, gB)))


def variant_single_5fold_C01():
    return ans(_rep(XA, yA, None, 5, repeats=1, C=0.1), _rep(XB, yB, gB, 5, repeats=1, C=0.1))


def variant_single_5fold_C10_noscale():
    return ans(_rep(XA, yA, None, 5, repeats=1, C=10, scale=False), _rep(XB, yB, gB, 5, repeats=1, C=10, scale=False))


def variant_A_half_B_honest():
    """Analyst reasons 'A is noise' and reports 0.5, computes B properly."""
    return ans(0.5, comp()["b_honest"])


TABLE = [
    ("correct (should PASS)", correct, True),
    ("variant 10-fold x3", variant_10fold, True),
    ("variant 3-fold x5", variant_3fold, True),
    ("variant leave-one-out/patient", variant_leave_one_out, True),
    ("variant single 5-fold C=0.1", variant_single_5fold_C01, True),
    ("variant single 5f C=10 no z", variant_single_5fold_C10_noscale, True),
    ("variant A=0.5 (noise), B honest", variant_A_half_B_honest, True),
    ("trap select-first (A)", trap_select_first_A, False),
    ("trap random-split (B)", trap_random_split_B, False),
    ("trap both", trap_both, False),
    ("trap chance 0.5/0.5", trap_chance, False),
    ("trap select-first grouped (B)", trap_select_first_B_grouped, False),
    ("trap resubstitution", trap_resubstitution, False),
]

if __name__ == "__main__":
    import json
    ev = HERE / "eval.json"
    truth = json.loads((HERE / "truth.json").read_text())
    print("constructed truth: A = 0.500 (exact), B = %.4f (Monte Carlo, 4000 fresh patients)" % truth["auc_dataset_b_fresh_patients"])
    print("components:", {k: round(v, 3) for k, v in comp().items()})
    ok = True
    for name, fn, expect in TABLE:
        a = fn()
        r = grade(ev, a)
        per = "".join("P" if g["passed"] else "f" for g in r["graders"])
        print(f"{name:34s} all={r['passed']!s:5s} [{per}] A={a['auc_dataset_a']:.2f} B={a['auc_dataset_b']:.2f}")
        ok &= (r["passed"] == expect)
    print("DISCRIMINATES (as designed)" if ok else "!! DOES NOT DISCRIMINATE")
    sys.exit(0 if ok else 1)
