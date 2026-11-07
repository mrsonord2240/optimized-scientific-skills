"""Synthetic task i1_ml_leakage: honest AUC of a filter-then-classify pipeline on omics-like data.

Everything here is FABRICATED.

Dataset A: 60 samples x 5000 features, features independent N(0,1), labels a random
           balanced assignment. The features carry no information about the label, so the
           expected AUC of ANY pipeline on new samples is exactly 0.5 (by construction).
Dataset B: 40 patients x 3 replicate samples x 2000 features. Label is per patient.
           - 8 features carry a real label effect (patient-level shift +-D_SIG),
           - 1200 features carry a large patient-specific offset shared by that patient's
             replicates (a "fingerprint", no relation to the label),
           - 792 features are pure replicate-level noise.
           The honest AUC on NEW patients is computed by Monte Carlo from the generating
           model: fit the specified pipeline on the released data, score 4000 freshly
           drawn patients. Independent of any CV code.

The graded pipeline: keep the 50 features with smallest two-sample t-test p-value
(sklearn f_classif, equivalent to pooled-variance t), z-score, L2 logistic regression
(sklearn default C=1).

Run:  python generate.py   -> data/*.csv, truth.json, eval.json   (deterministic)
"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
SEED_A, SEED_B = 17, 23                     # chosen by scan so common honest CV variants centre on the truth (NOTES.md)

P_A, N_A = 5000, 60
P_B, N_PAT, REPS = 2000, 40, 3
N_SIG, N_FP = 8, 1200
D_SIG, TAU_SIG, SIG_SIG = 1.1, 0.5, 1.0     # signal features: label shift, patient SD, replicate SD
TAU_FP, SIG_FP = 1.2, 0.4                   # fingerprint features: patient SD, replicate SD
SIG_NOISE = 1.0
N_FRESH = 4000                              # fresh patients for Monte Carlo truth

# Tolerance bands (absolute, asymmetric). See NOTES.md for rationale.
BAND_A = (0.20, 0.20)
BAND_B = (0.15, 0.13)


def make_pipeline_spec(C=1.0, k=50):
    return make_pipeline(SelectKBest(f_classif, k=k), StandardScaler(),
                         LogisticRegression(C=C, max_iter=5000))


def sim_a(rng):
    y = np.array([0, 1] * (N_A // 2))
    rng.shuffle(y)
    X = rng.normal(size=(N_A, P_A))
    return X, y


def draw_model(rng):
    """Column layout and signal signs (fixed for a 'world')."""
    perm = rng.permutation(P_B)
    sign = rng.choice([-1.0, 1.0], N_SIG)
    return dict(perm=perm, sign=sign)


def sim_b(rng, model, n_pat=N_PAT, reps=REPS, balanced=True):
    if balanced:
        y = np.array([0, 1] * (n_pat // 2))
        rng.shuffle(y)
    else:
        y = rng.integers(0, 2, n_pat)
    tau = np.zeros(P_B)
    sg = np.full(P_B, SIG_NOISE)
    beta = np.zeros(P_B)
    tau[:N_SIG], sg[:N_SIG], beta[:N_SIG] = TAU_SIG, SIG_SIG, model["sign"] * D_SIG
    tau[N_SIG:N_SIG + N_FP], sg[N_SIG:N_SIG + N_FP] = TAU_FP, SIG_FP
    pat = rng.normal(size=(n_pat, P_B)) * tau + np.outer(y - 0.5, beta)
    X = np.repeat(pat, reps, axis=0) + rng.normal(size=(n_pat * reps, P_B)) * sg
    X = X[:, np.argsort(model["perm"])]              # scatter blocks across columns
    return X, np.repeat(y, reps), np.repeat(np.arange(n_pat), reps)


def true_auc_b(model, X, y, seed):
    """AUC of the specified pipeline fitted on (X, y), on N_FRESH new patients (1 sample each)."""
    rng = np.random.default_rng(seed)
    Xf, yf, _ = sim_b(rng, model, n_pat=N_FRESH, reps=1, balanced=False)
    p = make_pipeline_spec().fit(X, y)
    return float(roc_auc_score(yf, p.decision_function(Xf)))


def cols(p):
    return [f"feat_{i + 1:04d}" for i in range(p)]


def build():
    """Return (dfA, dfB, truth dict, model). Pure function of SEED_A, SEED_B."""
    rng = np.random.default_rng([SEED_A, 0])
    # ---- A
    XA, yA = sim_a(rng)
    dfA = pd.DataFrame(XA, columns=cols(P_A))
    dfA.insert(0, "label", yA)
    dfA.insert(0, "sample_id", [f"A{i + 1:03d}" for i in range(N_A)])
    # ---- B (independent stream)
    rng = np.random.default_rng([SEED_B, 1])
    model = draw_model(rng)
    XB, yB, gB = sim_b(rng, model)
    order = rng.permutation(len(yB))                 # shuffle rows: replicates not adjacent
    XB, yB, gB = XB[order], yB[order], gB[order]
    pid = np.array([f"PT{g + 1:02d}" for g in gB])
    dfB = pd.DataFrame(XB, columns=cols(P_B))
    dfB.insert(0, "label", yB)
    dfB.insert(0, "patient_id", pid)
    dfB.insert(0, "sample_id", [f"B{i + 1:03d}" for i in range(len(yB))])
    # ---- truth, computed on exactly what is written to disk (4-dp rounding)
    XB_disk = np.round(XB, 4)
    auc_b = true_auc_b(model, XB_disk, yB, SEED_B + 1000)
    signal_cols = sorted(int(np.where(model["perm"] == j)[0][0]) + 1 for j in range(N_SIG))
    truth = dict(
        auc_dataset_a=0.5,
        auc_dataset_b_fresh_patients=round(auc_b, 4),
        note="A: labels independent of features (exact). B: Monte Carlo on 4000 fresh patients from the generating model.",
        signal_feature_numbers_B=signal_cols,
        n_fingerprint_features_B=N_FP,
        D_SIG=D_SIG, TAU_FP=TAU_FP,
    )
    return dfA, dfB, truth, model


TASK = (
    "A collaborator has sent two small omics datasets, each with a binary phenotype (label: 1 = case, 0 = control). "
    "Values are normalised, log-scale measurements.\n\n"
    "- data/dataset_a.csv: 60 samples (one row each). Columns: sample_id, label, then 5,000 numeric features (feat_0001 ... feat_5000).\n"
    "- data/dataset_b.csv: 120 samples from 40 patients, 3 samples per patient. Columns: sample_id, patient_id, label "
    "(the label belongs to the patient), then 2,000 numeric features (feat_0001 ... feat_2000).\n\n"
    "For each dataset we intend to build the same classifier: keep the 50 features with the smallest p-values from a "
    "two-sample t-test between the two classes, z-score them, and fit an L2-regularised logistic regression "
    "(scikit-learn LogisticRegression, default C = 1) on those 50 features.\n\n"
    "We need a trustworthy estimate of how well this pipeline will perform in practice. For each dataset, report the "
    "expected ROC AUC of the whole pipeline (feature selection included) on samples (dataset A) or patients (dataset B) "
    "that were not used for any fitting or selection step. Report each AUC to 2 decimal places.\n\n"
    "Return EXACTLY:\n\n<EVAL_ANSWER>\n{\"auc_dataset_a\": 0.0, \"auc_dataset_b\": 0.0}\n</EVAL_ANSWER>"
)


def graders(truth):
    ta, tb = truth["auc_dataset_a"], round(truth["auc_dataset_b_fresh_patients"], 2)
    def one(field, t, band):
        return {"type": "numeric_tolerance", "config": {
            "tolerances": {field: {"type": "absolute", "lower": band[0], "upper": band[1]}},
            "ground_truth": {field: t}}}
    return [one("auc_dataset_a", ta, BAND_A), one("auc_dataset_b", tb, BAND_B)]


def main():
    warnings.filterwarnings("ignore")
    (HERE / "data").mkdir(exist_ok=True)
    dfA, dfB, truth, _ = build()
    dfA.to_csv(HERE / "data" / "dataset_a.csv", index=False, float_format="%.4f", lineterminator="\n")
    dfB.to_csv(HERE / "data" / "dataset_b.csv", index=False, float_format="%.4f", lineterminator="\n")
    (HERE / "truth.json").write_text(json.dumps(truth, indent=2) + "\n", encoding="utf-8")
    ev = {"id": "i1_ml_leakage", "replicates": None, "synthetic": True, "task": TASK,
          "graders": graders(truth), "data_dir": "data",
          "relevant_skills": ["bio-machine-learning-model-validation"],
          "design_source": "Feature selection outside CV (Ambroise & McLachlan 2002, PNAS 99:6562) and replicate samples from one patient split across folds (subject-wise vs record-wise CV)."}
    (HERE / "eval.json").write_text(json.dumps(ev, indent=2) + "\n", encoding="utf-8")
    print("truth:", truth["auc_dataset_a"], truth["auc_dataset_b_fresh_patients"])


if __name__ == "__main__":
    main()
