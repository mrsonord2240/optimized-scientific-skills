"""Reference solvers for the synthetic microbiome differential-abundance task.

`correct()`   : compositionally aware analysis (CLR + linear model, adjusted for sex, BH 5%).
`variant_*()` : other defensible compositional methods (ALDEx2-style, median-centred, LinDA-like);
                must also PASS (they show the tolerances admit legitimate method variants).
`trap_*()`    : the naive relative-abundance analyses. They must FAIL.

Usage:
  python solvers.py            run every solver through the frozen graders (exit 0 if as designed)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from grading import grade  # noqa: E402

D = HERE / "data"
counts = pd.read_csv(D / "counts.tsv", sep="\t", index_col=0)
meta = pd.read_csv(D / "metadata.tsv", sep="\t").set_index("sample").loc[counts.columns]
X = counts.values.astype(float)                 # taxa x samples
ids = np.array(counts.index)
case = (meta["group"] == "case").values.astype(float)
female = (meta["sex"] == "F").values.astype(float)
depth = X.sum(0)
MIN_PREV = 8                                    # keep taxa detected in >= 8 of 40 samples
ALPHA = 0.05


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def ols_group(Y, covar_female=True):
    """Row-wise OLS of Y (taxa x samples) on [1, case, (female)]. Returns beta_case, se, p, df."""
    cols = [np.ones_like(case), case] + ([female] if covar_female else [])
    M = np.column_stack(cols)
    XtXi = np.linalg.inv(M.T @ M)
    B = Y @ M @ XtXi                                        # taxa x p
    R = Y - B @ M.T
    df = M.shape[0] - M.shape[1]
    s2 = (R ** 2).sum(1) / df
    se = np.sqrt(s2 * XtXi[1, 1])
    with np.errstate(divide="ignore", invalid="ignore"):
        t = B[:, 1] / se
    p = 2 * stats.t.sf(np.abs(t), df)
    p = np.where(np.isfinite(p), p, 1.0)
    return B[:, 1], se, p, df


def keep_mask(minprev=MIN_PREV):
    return (X > 0).sum(1) >= minprev


def answer(sig):
    return {"differential_taxa": sig}


# ---------------------------------------------------------------- compositional methods
def clr(Xm, pc=0.5):
    L = np.log(Xm + pc)
    return L - L.mean(0, keepdims=True)


def _clr_ols(pc=0.5, adjust=True, minprev=MIN_PREV):
    keep = keep_mask(minprev)
    Y = clr(X[keep], pc)                        # CLR over the retained taxa
    _, _, p, _ = ols_group(Y, adjust)
    q = bh(p)
    return sorted(ids[keep][q < ALPHA].tolist())


def correct():
    """CLR (pseudocount 0.5) + OLS on group and sex, prevalence filter, BH 5%."""
    return answer(_clr_ols(0.5, True))


def variant_clr_unadjusted():
    return answer(_clr_ols(0.5, False))


def variant_clr_pc1_noprev():
    return answer(_clr_ols(1.0, True, minprev=2))


def variant_median_centred():
    """Reference-frame centring: log(x+1) minus the per-sample median log-ratio of detected taxa."""
    keep = keep_mask()
    L = np.log(X[keep] + 1.0)
    det = X[keep] > 0
    med = np.array([np.median(L[det[:, j], j]) for j in range(L.shape[1])])
    _, _, p, _ = ols_group(L - med, True)
    return answer(sorted(ids[keep][bh(p) < ALPHA].tolist()))


def variant_aldex_style(nmc=128, seed=1):
    """ALDEx2-like: Dirichlet Monte-Carlo instances of (count + 0.5), CLR, per-instance test, mean p, BH."""
    rng = np.random.default_rng(seed)
    keep = keep_mask()
    Xk = X[keep]
    ps = np.zeros((nmc, Xk.shape[0]))
    for m in range(nmc):
        g = rng.gamma(Xk + 0.5, 1.0)
        pr = g / g.sum(0, keepdims=True)
        L = np.log(pr)
        Y = L - L.mean(0, keepdims=True)
        ps[m] = ols_group(Y, True)[2]
    return answer(sorted(ids[keep][bh(ps.mean(0)) < ALPHA].tolist()))


def variant_linda_like():
    """LinDA-style: OLS on log(count+0.5) (adjusted), then remove the mode of the coefficients (compositional bias)."""
    keep = keep_mask()
    Y = np.log(X[keep] + 0.5)
    b, se, _, df = ols_group(Y, True)
    kde = stats.gaussian_kde(b)
    grid = np.linspace(b.min(), b.max(), 2000)
    mode = grid[np.argmax(kde(grid))]
    t = (b - mode) / se
    p = 2 * stats.t.sf(np.abs(t), df)
    return answer(sorted(ids[keep][bh(p) < ALPHA].tolist()))


def variant_deseq_style_nb():
    """Median-of-ratios size factors (DESeq2/'poscounts'-like) + log1p-normalised OLS."""
    keep = keep_mask()
    Xk = X[keep] + 0.5
    lg = np.log(Xk)
    ref = lg.mean(1, keepdims=True)
    sf = np.exp(np.median(lg - ref, axis=0))
    Y = np.log(Xk / sf)
    return answer(sorted(ids[keep][bh(ols_group(Y, True)[2]) < ALPHA].tolist()))


# ---------------------------------------------------------------- naive relative-abundance traps
def _prop():
    return X / depth


def trap_wilcoxon_relabund():
    keep = keep_mask()
    P = _prop()[keep]
    p = np.array([stats.mannwhitneyu(r[case == 1], r[case == 0], alternative="two-sided").pvalue for r in P])
    return answer(sorted(ids[keep][bh(p) < ALPHA].tolist()))


def trap_ttest_relabund():
    keep = keep_mask()
    P = _prop()[keep]
    p = np.array([stats.ttest_ind(r[case == 1], r[case == 0], equal_var=False).pvalue for r in P])
    p = np.where(np.isfinite(p), p, 1.0)
    return answer(sorted(ids[keep][bh(p) < ALPHA].tolist()))


def trap_rarefy_wilcoxon(seed=7):
    rng = np.random.default_rng(seed)
    dmin = int(depth.min())
    R = np.zeros_like(X)
    for j in range(X.shape[1]):
        R[:, j] = rng.multivariate_hypergeometric(X[:, j].astype(int), dmin)
    keep = (R > 0).sum(1) >= MIN_PREV
    p = np.array([stats.mannwhitneyu(r[case == 1], r[case == 0], alternative="two-sided").pvalue for r in R[keep]])
    return answer(sorted(ids[keep][bh(p) < ALPHA].tolist()))


def trap_logcpm_ols():
    """Library-size (TSS/CPM) normalisation, log transform, linear model adjusted for sex (limma-like)."""
    keep = keep_mask()
    Y = np.log2(_prop()[keep] * 1e6 + 1.0)
    return answer(sorted(ids[keep][bh(ols_group(Y, True)[2]) < ALPHA].tolist()))


def trap_nb_glm_offset():
    """Negative-binomial GLM with log(library size) offset, group + sex (edgeR/DESeq2-with-TSS-like)."""
    import statsmodels.api as sm
    keep = keep_mask()
    M = np.column_stack([np.ones(len(case)), case, female])
    ps = []
    for r in X[keep]:
        try:
            # moment estimate of NB dispersion from a quasi-Poisson fit, then NB GLM
            pf = sm.GLM(r, M, family=sm.families.Poisson(), offset=np.log(depth)).fit()
            mu = pf.mu
            alpha = max(((r - mu) ** 2 - mu).sum() / (mu ** 2).sum(), 1e-3)
            nb = sm.GLM(r, M, family=sm.families.NegativeBinomial(alpha=alpha), offset=np.log(depth)).fit()
            ps.append(nb.pvalues[1])
        except Exception:
            ps.append(1.0)
    return answer(sorted(ids[keep][bh(np.nan_to_num(np.array(ps), nan=1.0)) < ALPHA].tolist()))


def variant_clr_wilcoxon():
    keep = keep_mask()
    Y = clr(X[keep], 0.5)
    p = np.array([stats.mannwhitneyu(r[case == 1], r[case == 0], alternative="two-sided").pvalue for r in Y])
    return answer(sorted(ids[keep][bh(p) < ALPHA].tolist()))


def trap_wilcoxon_relabund_nofilter():
    """Wilcoxon on proportions with every non-degenerate taxon tested (all-zero rows get p=1)."""
    P = _prop()
    p = np.ones(len(P))
    for i, r in enumerate(P):
        if r.max() > r.min():
            p[i] = stats.mannwhitneyu(r[case == 1], r[case == 0], alternative="two-sided").pvalue
    return answer(sorted(ids[bh(p) < ALPHA].tolist()))


SOLVERS = [  # (label, fn, expected_to_pass)
    ("correct: CLR+OLS(sex)", correct, True),
    ("variant CLR unadjusted", variant_clr_unadjusted, True),
    ("variant CLR pc=1 no-filter", variant_clr_pc1_noprev, True),
    ("variant median-centred", variant_median_centred, True),
    ("variant ALDEx2-style", variant_aldex_style, True),
    ("variant LinDA-like", variant_linda_like, True),
    ("variant CLR + Wilcoxon", variant_clr_wilcoxon, True),
    ("variant DESeq-sf (median ratio)", variant_deseq_style_nb, True),
    ("trap Wilcoxon on proportions", trap_wilcoxon_relabund, False),
    ("trap Wilcoxon, no filter", trap_wilcoxon_relabund_nofilter, False),
    ("trap Welch t on proportions", trap_ttest_relabund, False),
    ("trap rarefy + Wilcoxon", trap_rarefy_wilcoxon, False),
    ("trap log-CPM OLS(sex)", trap_logcpm_ols, False),
    ("trap NB-GLM lib-size offset", trap_nb_glm_offset, False),
]


if __name__ == "__main__":
    truth = set(json.loads((HERE / "truth.json").read_text())["truth_taxa"])
    ev = HERE / "eval.json"
    ok = True
    print(f"{'solver':34s} n  TP FP FN  prec  rec  graders all")
    for name, fn, expect in SOLVERS:
        a = fn()
        s = set(a["differential_taxa"])
        tp, fp, fn_ = len(s & truth), len(s - truth), len(truth - s)
        r = grade(ev, a)
        per = "".join("P" if g["passed"] else "f" for g in r["graders"])
        prec = tp / len(s) if s else 0.0
        print(f"{name:34s} {len(s):3d} {tp:2d} {fp:3d} {fn_:2d} {prec:5.2f} {tp / len(truth):4.2f}  [{per}] {r['passed']}")
        ok &= (r["passed"] == expect)
    print("DISCRIMINATES (as designed)" if ok else "!! DOES NOT DISCRIMINATE")
    sys.exit(0 if ok else 1)
