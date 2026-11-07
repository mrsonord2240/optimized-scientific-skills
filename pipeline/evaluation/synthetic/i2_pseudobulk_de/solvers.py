"""Reference solvers for i2_pseudobulk_de.

`correct*`  : pseudobulk by donor x condition after excluding donors with < 5 cells in either condition,
              donor as blocking factor (paired). Three variants: OLS on log2 CPM, quasi-Poisson GLM, paired t-test.
`trap_*`    : cell-level Wilcoxon, cell-level Welch t-test, pseudobulk ignoring donor pairing,
              pseudobulk without the donor-size filter.
Usage: python solvers.py
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from grading import grade  # noqa: E402

D = HERE / "data"
MIN_CELLS, MIN_TOTAL, ALPHA, LFC = 5, 20, 0.05, 1.0


def load():
    with gzip.open(D / "matrix.mtx.gz", "rt") as f:
        M = sp.coo_matrix(  # genes x cells
            _read_mtx(f)).tocsr()
    feats = pd.read_csv(D / "features.tsv.gz", sep="\t", header=None, names=["id", "sym", "type"])
    bc = pd.read_csv(D / "barcodes.tsv.gz", header=None)[0].tolist()
    cells = pd.read_csv(D / "cells.tsv", sep="\t").set_index("barcode").loc[bc].reset_index()
    return M.T.tocsr(), feats.sym.values, cells      # cells x genes


def _read_mtx(f):
    from scipy.io import mmread
    return mmread(f)


X, GENES, CELLS = load()


def bh(p):
    p = np.asarray(p, float)
    out = np.ones_like(p)
    ok = ~np.isnan(p)
    q = p[ok]
    o = np.argsort(q)
    r = q[o] * len(q) / (np.arange(len(q)) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    res = np.empty_like(q)
    res[o] = np.minimum(r, 1)
    out[ok] = res
    return out


def retained_donors(cells=CELLS):
    t = pd.crosstab(cells.donor, cells.condition)
    return sorted(t.index[(t.control >= MIN_CELLS) & (t.treated >= MIN_CELLS)])


def pick(res_p, res_lfc, mask):
    padj = np.full(len(res_p), np.nan)
    padj[mask] = bh(res_p[mask])
    sel = mask & (padj < ALPHA) & (np.abs(res_lfc) > LFC)
    return [g for g, s in zip(GENES, sel) if s]


# ----- pseudobulk ------------------------------------------------------------------------
def pseudobulk(donors):
    c = CELLS[CELLS.donor.isin(donors)]
    keys = c.groupby(["donor", "condition"]).indices
    names = sorted(keys)
    Xr = X[c.index.values]
    PB = np.vstack([np.asarray(Xr[keys[k]].sum(0)).ravel() for k in names])
    meta = pd.DataFrame(names, columns=["donor", "condition"])
    return PB, meta


def gene_mask(donors):
    keep = CELLS.donor.isin(donors).values
    return np.asarray(X[keep].sum(0)).ravel() >= MIN_TOTAL


def pb_logcpm(PB, prior=1.0):
    return np.log2(PB / PB.sum(1, keepdims=True) * 1e6 + prior)


def ols_pb(donors, paired=True, prior=1.0):
    PB, meta = pseudobulk(donors)
    Y = pb_logcpm(PB, prior)
    tr = (meta.condition == "treated").astype(float).values
    cols = [np.ones(len(meta)), tr]
    if paired:
        for d in sorted(set(meta.donor))[1:]:
            cols.append((meta.donor == d).astype(float).values)
    A = np.column_stack(cols)
    beta, *_ = np.linalg.lstsq(A, Y, rcond=None)
    res = Y - A @ beta
    df = len(meta) - A.shape[1]
    s2 = (res ** 2).sum(0) / df
    cov = np.linalg.inv(A.T @ A)[1, 1]
    t = beta[1] / np.sqrt(s2 * cov + 1e-30)
    p = 2 * stats.t.sf(np.abs(t), df)
    return p, beta[1], gene_mask(donors)


def paired_t_pb(donors, prior=1.0):
    PB, meta = pseudobulk(donors)
    Y = pb_logcpm(PB, prior)
    ds = sorted(set(meta.donor))
    ic = [meta.index[(meta.donor == d) & (meta.condition == "control")][0] for d in ds]
    it = [meta.index[(meta.donor == d) & (meta.condition == "treated")][0] for d in ds]
    diff = Y[it] - Y[ic]
    t, p = stats.ttest_1samp(diff, 0.0, axis=0)
    return p, diff.mean(0), gene_mask(donors)


def qpoisson_pb(donors, paired=True):
    import statsmodels.api as sm
    PB, meta = pseudobulk(donors)
    lib = PB.sum(1)
    tr = (meta.condition == "treated").astype(float).values
    cols = [np.ones(len(meta)), tr]
    if paired:
        for d in sorted(set(meta.donor))[1:]:
            cols.append((meta.donor == d).astype(float).values)
    A = np.column_stack(cols)
    mask = gene_mask(donors)
    p = np.full(len(GENES), np.nan)
    lf = np.zeros(len(GENES))
    off = np.log(lib / np.exp(np.log(lib).mean()))
    import warnings
    warnings.filterwarnings("ignore")
    for j in np.where(mask)[0]:
        y = PB[:, j]
        try:
            r = sm.GLM(y, A, family=sm.families.Poisson(), offset=off).fit(scale="X2")
            lf[j] = r.params[1] / np.log(2)
            p[j] = r.pvalues[1]
        except Exception:
            p[j] = 1.0
    return p, lf, mask


# ----- cell level ------------------------------------------------------------------------
def lognorm(sub):
    lib = np.asarray(sub.sum(1)).ravel()
    Xn = sp.diags(1e4 / lib) @ sub
    return Xn


def cell_level(donors, test):
    keep = CELLS.donor.isin(donors).values if donors is not None else np.ones(len(CELLS), bool)
    sub = X[keep]
    cond = (CELLS.condition.values[keep] == "treated")
    Xn = lognorm(sub).tocsc()
    mask = np.asarray(sub.sum(0)).ravel() >= MIN_TOTAL
    L = Xn.copy()
    L.data = np.log1p(L.data)
    Ld = L.toarray()
    t, c = Ld[cond], Ld[~cond]
    if test == "wilcoxon":
        r = stats.rankdata(Ld, axis=0)
        n1, n2 = cond.sum(), (~cond).sum()
        N = n1 + n2
        U = r[cond].sum(0) - n1 * (n1 + 1) / 2
        mu = n1 * n2 / 2
        # tie-corrected variance
        tie = np.zeros(Ld.shape[1])
        for j in np.where(mask)[0]:
            _, cnt = np.unique(Ld[:, j], return_counts=True)
            tie[j] = (cnt ** 3 - cnt).sum()
        sd = np.sqrt(n1 * n2 / 12 * ((N + 1) - tie / (N * (N - 1))) + 1e-30)
        z = (U - mu) / sd
        p = 2 * stats.norm.sf(np.abs(z))
    else:
        tt, p = stats.ttest_ind(t, c, axis=0, equal_var=False)
    # scanpy-style log2FC from mean of expm1 of log-normalised values
    mt, mc = np.expm1(t).mean(0), np.expm1(c).mean(0)
    lfc = np.log2((mt + 1e-9) / (mc + 1e-9))
    return np.nan_to_num(p, nan=1.0), lfc, mask


# ----- solvers -----------------------------------------------------------------------------
def _ans(genes):
    return {"de_genes": genes}


def correct():                                   # primary reference: OLS on log2CPM, donor as covariate
    return _ans(pick(*ols_pb(retained_donors())))


def correct_qp():
    return _ans(pick(*qpoisson_pb(retained_donors())))


def correct_paired_t():
    return _ans(pick(*paired_t_pb(retained_donors())))


def trap_wilcoxon_cells():
    return _ans(pick(*cell_level(retained_donors(), "wilcoxon")))


def trap_wilcoxon_cells_alldonors():
    return _ans(pick(*cell_level(None, "wilcoxon")))


def trap_ttest_cells():
    return _ans(pick(*cell_level(retained_donors(), "ttest")))


def trap_pseudobulk_unpaired():
    return _ans(pick(*ols_pb(retained_donors(), paired=False)))


def trap_pseudobulk_unpaired_qp():
    return _ans(pick(*qpoisson_pb(retained_donors(), paired=False)))


def trap_pseudobulk_no_filter():
    alld = sorted(set(CELLS.donor))
    # donors with a missing condition would break pairing; all donors here have cells in both conditions
    return _ans(pick(*ols_pb(alld)))


def trap_pseudobulk_no_filter_qp():
    alld = sorted(set(CELLS.donor))
    return _ans(pick(*qpoisson_pb(alld)))


def trap_pseudobulk_no_filter_paired_t():
    alld = sorted(set(CELLS.donor))
    return _ans(pick(*paired_t_pb(alld)))


# (name, fn, expected to pass?)
SOLVERS = [
    ("correct OLS log2CPM+donor", correct, True),
    ("correct quasi-Poisson+donor", correct_qp, True),
    ("correct paired t-test", correct_paired_t, True),
    ("trap cell Wilcoxon", trap_wilcoxon_cells, False),
    ("trap cell Wilcoxon (all donors)", trap_wilcoxon_cells_alldonors, False),
    ("trap cell Welch t", trap_ttest_cells, False),
    ("trap pseudobulk unpaired OLS", trap_pseudobulk_unpaired, False),
    ("trap pseudobulk unpaired QP", trap_pseudobulk_unpaired_qp, False),
    ("trap pseudobulk no filter OLS", trap_pseudobulk_no_filter, False),
    # BLIND SPOT: count-based GLM weights tiny samples by their Poisson variance, so the filter barely matters
    ("BLIND SPOT no filter QP-GLM", trap_pseudobulk_no_filter_qp, True),
    ("trap pseudobulk no filter t", trap_pseudobulk_no_filter_paired_t, False),
]


if __name__ == "__main__":
    import json
    truth = set(json.loads((HERE / "truth.json").read_text())["de_genes"])
    ev = HERE / "eval.json"
    only = sys.argv[1:]
    ok = True
    print(f"{'solver':34s} {'n':>4s} {'TP':>4s} {'prec':>6s} {'rec':>6s}  grade  expect")
    for name, fn, expect in SOLVERS:
        if only and not any(o in name for o in only):
            continue
        a = fn()
        g = set(a["de_genes"])
        tp = len(g & truth)
        r = grade(ev, a)
        per = "".join("P" if x["passed"] else "f" for x in r["graders"])
        prec = tp / max(len(a["de_genes"]), 1)
        print(f"{name:34s} {len(g):4d} {tp:4d} {prec:6.3f} {tp / len(truth):6.3f}  [{per}]  {'PASS' if expect else 'FAIL'}"
              f"{'' if r['passed'] == expect else '   <-- MISMATCH'}", flush=True)
        ok &= (r["passed"] == expect)
    print("DISCRIMINATES (as designed)" if ok else "!! DOES NOT DISCRIMINATE")
    sys.exit(0 if ok else 1)
