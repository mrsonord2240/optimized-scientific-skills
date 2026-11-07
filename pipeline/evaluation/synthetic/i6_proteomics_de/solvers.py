"""Reference solvers for the synthetic proteomics differential-abundance task.

`correct*` : log2, per-sample median normalisation, eligibility filter, per-protein linear model
             (condition + batch) on observed values only (no imputation), moderated (limma-style
             empirical Bayes) or ordinary t, BH at 5%. `correct_perseus` adds MNAR-aware
             (downshifted-Gaussian) imputation before a moderated test.
`trap_*`   : the shortcuts named in the brief. They must FAIL.

Usage:
  python solvers.py                # run everything through the graders in eval.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from grading import grade  # noqa: E402

D = HERE / "data"
MIN_VALID = 3
ALPHA = 0.05

tab = pd.read_csv(D / "protein_intensities.tsv", sep="\t", dtype={"Reverse": str, "Potential contaminant": str},
                  keep_default_na=False, na_values=[""])
smp = pd.read_csv(D / "samples.tsv", sep="\t")
tab = tab[tab["Reverse"].isna() & tab["Potential contaminant"].isna()].reset_index(drop=True)
IDS = tab["Protein IDs"].tolist()
RAW = tab[[f"Intensity {s}" for s in smp["sample"]]].to_numpy(float)     # NaN = missing
TREAT = (smp["condition"] == "treated").to_numpy()
B2 = (smp["batch"] == 2).to_numpy()


# ---------------------------------------------------------------- helpers
def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def eligible(mat):
    v = ~np.isnan(mat)
    return (v[:, TREAT].sum(1) >= MIN_VALID) & (v[:, ~TREAT].sum(1) >= MIN_VALID)


def log2n(raw=RAW, normalise=True):
    y = np.log2(np.where(raw > 0, raw, np.nan))
    if normalise:
        y = y - np.nanmedian(y, axis=0)[None, :]
    return y


def trigamma_inv(x):
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = polygamma(1, y)
        dy = tri * (1 - tri / x) / polygamma(2, y)
        y = y + dy
        if abs(dy) < 1e-10:
            break
    return y


def fit_fdist(s2, df):
    """limma fitFDist (method of moments) for prior d0, s0^2; df may vary."""
    z = np.log(s2)
    e = z - digamma(df / 2) + np.log(df / 2)
    emean = e.mean()
    evar = e.var(ddof=1) - np.mean(polygamma(1, df / 2))
    if evar > 0:
        d0 = 2 * trigamma_inv(evar)
        s02 = np.exp(emean + digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0, s02 = np.inf, np.exp(emean)
    return d0, s02


def lm_test(y, use_batch=True, moderate=True, elig=None):
    """Per-protein OLS of y on [1, treated, batch2]; observed values only. Returns p-values (NaN if not testable)."""
    n, m = y.shape
    X = np.column_stack([np.ones(m), TREAT.astype(float)] + ([B2.astype(float)] if use_batch else []))
    beta = np.full(n, np.nan)
    s2 = np.full(n, np.nan)
    dfr = np.full(n, np.nan)
    v = np.full(n, np.nan)
    for i in range(n):
        if elig is not None and not elig[i]:
            continue
        ok = ~np.isnan(y[i])
        Xi, yi = X[ok], y[i, ok]
        if np.linalg.matrix_rank(Xi) < X.shape[1]:
            continue
        XtXi = np.linalg.inv(Xi.T @ Xi)
        b = XtXi @ Xi.T @ yi
        r = yi - Xi @ b
        d = ok.sum() - X.shape[1]
        if d < 1:
            continue
        beta[i], s2[i], dfr[i], v[i] = b[1], (r @ r) / d, d, XtXi[1, 1]
    t_ok = ~np.isnan(beta)
    if moderate:
        d0, s02 = fit_fdist(np.maximum(s2[t_ok], 1e-8), dfr[t_ok])
        s2m = (d0 * s02 + dfr * s2) / (d0 + dfr)
        dft = dfr + d0
    else:
        s2m, dft = s2, dfr
    t = beta / np.sqrt(s2m * v)
    p = 2 * stats.t.sf(np.abs(t), dft)
    return beta, p


def call(mat, p):
    q = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    q[ok] = bh(p[ok])
    return [IDS[i] for i in np.where(q <= ALPHA)[0]]


# ---------------------------------------------------------------- correct variants
def correct():
    y = log2n()
    el = eligible(y)
    _, p = lm_test(y, use_batch=True, moderate=True, elig=el)
    return {"differentially_abundant_proteins": call(y, p)}


def correct_ols():
    y = log2n()
    el = eligible(y)
    _, p = lm_test(y, use_batch=True, moderate=False, elig=el)
    return {"differentially_abundant_proteins": call(y, p)}


def correct_perseus(seed=0):
    """Left-censored (MNAR) imputation: downshifted Gaussian per sample (Perseus defaults 1.8 SD / 0.3 width)."""
    rng = np.random.default_rng(seed)
    y = log2n()
    el = eligible(y)
    yi = y.copy()
    for j in range(y.shape[1]):
        c = y[:, j]
        mu, sd = np.nanmean(c), np.nanstd(c)
        m = np.isnan(c)
        yi[m, j] = rng.normal(mu - 1.8 * sd, 0.3 * sd, m.sum())
    _, p = lm_test(yi, use_batch=True, moderate=True, elig=el)
    return {"differentially_abundant_proteins": call(y, p)}


# ---------------------------------------------------------------- traps
def _student(y, elig, welch=False):
    p = np.full(y.shape[0], np.nan)
    for i in np.where(elig)[0]:
        a, b = y[i, TREAT], y[i, ~TREAT]
        a, b = a[~np.isnan(a)], b[~np.isnan(b)]
        if len(a) > 1 and len(b) > 1:
            p[i] = stats.ttest_ind(a, b, equal_var=not welch).pvalue
    return p


def trap_no_norm_min_impute(per_protein=True):
    """No normalisation; missing -> minimum observed value (per protein, or global); Student t; BH."""
    y = log2n(normalise=False)
    el = eligible(y)
    yi = y.copy()
    if per_protein:
        mn = np.nanmin(np.where(el[:, None], y, 0.0), axis=1)
        for j in range(y.shape[1]):
            yi[np.isnan(yi[:, j]), j] = mn[np.isnan(yi[:, j])]
    else:
        yi[np.isnan(yi)] = np.nanmin(y)
    p = _student(yi, el)
    return {"differentially_abundant_proteins": call(y, p)}


def trap_zero_impute(normalise=True):
    """Missing -> 0 intensity, log2(x+1), (per-sample median normalisation from observed values), Student t; BH."""
    z = np.log2(np.where(np.isnan(RAW), 0.0, RAW) + 1.0)
    if normalise:
        med = np.nanmedian(log2n(normalise=False), axis=0)
        z = np.where(np.isnan(RAW), 0.0, z - med[None, :])
    el = eligible(RAW)
    p = _student(z, el)
    return {"differentially_abundant_proteins": call(z, p)}


def trap_min_impute_normalised():
    """Isolates the imputation: median-normalised, per-protein-min imputation, Student t; BH."""
    y = log2n()
    el = eligible(y)
    mn = np.nanmin(np.where(el[:, None], y, 0.0), axis=1)
    yi = np.where(np.isnan(y), mn[:, None], y)
    return {"differentially_abundant_proteins": call(y, _student(yi, el))}


def extra_student_norm_observed():
    """Median-normalised, observed values only, Student t (no batch, no moderation), BH."""
    y = log2n()
    return {"differentially_abundant_proteins": call(y, _student(y, eligible(y)))}


def trap_ignore_batch():
    """Median-normalised, observed values only, moderated test with condition only (batch ignored)."""
    y = log2n()
    el = eligible(y)
    _, p = lm_test(y, use_batch=False, moderate=True, elig=el)
    return {"differentially_abundant_proteins": call(y, p)}


def trap_welch_complete_no_correction():
    """Median-normalised; complete cases only (all 12 observed); Welch t; raw p < 0.05, no multiplicity correction."""
    y = log2n()
    cc = ~np.isnan(y).any(axis=1)
    p = _student(y, cc, welch=True)
    return {"differentially_abundant_proteins": [IDS[i] for i in np.where(p < ALPHA)[0]]}


# ---------------------------------------------------------------- run
def _metrics(ans, truth):
    got, tr = {x.lower() for x in ans["differentially_abundant_proteins"]}, {x.lower() for x in truth}
    tp = len(got & tr)
    return len(got), tp, (tp / len(got) if got else 0.0), tp / len(tr)


if __name__ == "__main__":
    truth = json.loads((HERE / "truth.json").read_text())["de_eligible"]
    THRESH_CHECK = {"precision_at_k": 0.85, "recall_at_k": 0.80}
    ev = HERE / "eval.json"
    spec = json.loads(ev.read_text(encoding="utf-8"))
    assert spec["graders"][0]["config"]["scoring"]["pass_thresholds"] == THRESH_CHECK
    rows = [
        ("correct (mod. t + batch)", correct, True),
        ("correct (OLS t + batch)", correct_ols, True),
        ("info: Perseus-style impute+mod", correct_perseus, "info"),
        ("trap no-norm+min imp+t", trap_no_norm_min_impute, False),
        ("  (global-min variant)", lambda: trap_no_norm_min_impute(False), False),
        ("info: norm + min imp + t", trap_min_impute_normalised, "info"),
        ("trap zero imputation", trap_zero_impute, False),
        ("  (zero imp., no norm)", lambda: trap_zero_impute(False), False),
        ("trap ignore batch (mod. t)", trap_ignore_batch, "blind"),  # named trap: outcome recorded, not asserted
        ("trap Welch cc no correction", trap_welch_complete_no_correction, False),
        ("info: norm+Student, no batch", extra_student_norm_observed, "info"),
    ]
    ok = True
    print(f"truth (eligible true-DE): {len(truth)}")
    print(f"{'solver':29s} {'pass':5s} n_rep  TP  prec   rec   grader")
    for name, fn, expect in rows:
        a = fn()
        r = grade(ev, a)
        n, tp, pr, rc = _metrics(a, truth)
        per = "".join("P" if g["passed"] else "f" for g in r["graders"])
        tag = ""
        if expect == "blind":
            tag = "  <-- BLIND SPOT (named trap passes)" if r["passed"] else "  (fails)"
        elif expect == "info":
            tag = "  (informational)"
        elif r["passed"] != expect:
            ok = False
            tag = "  !! unexpected"
        print(f"{name:29s} {r['passed']!s:5s} {n:5d} {tp:4d} {pr:5.3f} {rc:5.3f}  [{per}]{tag}")
    print("DISCRIMINATES (as designed)" if ok else "!! DOES NOT DISCRIMINATE")
    sys.exit(0 if ok else 1)
