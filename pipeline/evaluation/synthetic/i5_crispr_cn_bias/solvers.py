"""Reference solvers for the synthetic CRISPR copy-number-bias task.

`correct*`  : copy-number-aware calls (three legitimate variants) -> must PASS.
`trap_*`    : shortcuts that ignore copy number -> must FAIL.
All calls: gene-level (median of 4 guides) LFC <= -1 and BH-FDR < 0.05 vs the
non-targeting-control guide distribution (task convention), unless stated.

Usage: python solvers.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from grading import grade  # noqa: E402

D = HERE / "data"
counts = pd.read_csv(D / "counts.tsv", sep="\t", index_col=0)
lib = pd.read_csv(D / "library.tsv", sep="\t", index_col=0)
coords = pd.read_csv(D / "gene_coordinates.tsv", sep="\t")
segs = pd.read_csv(D / "copy_number_segments.tsv", sep="\t")

FLOOR = -1.0
ALPHA = 0.05


# ---- basics ----------------------------------------------------------------
def guide_lfc():
    cpm = counts / counts.sum() * 1e6
    lfc = np.log2((cpm[["late_rep1", "late_rep2"]] + 0.5).div(cpm["plasmid"] + 0.5, axis=0)).mean(axis=1)
    lfc = lfc - lfc[lib.gene == "NTC"].median()
    return lfc.rename("lfc")


def gene_cn():
    out = {}
    for _, r in coords.iterrows():
        mid = (r.start + r.end) / 2
        s = segs[(segs.chromosome == r.chromosome) & (segs.start <= mid) & (segs.end >= mid)]
        out[r.gene] = float(s.copy_number.iloc[0])
    return pd.Series(out, name="cn")


def bh(p):
    p = np.asarray(p)
    o = np.argsort(p)
    q = np.empty_like(p)
    ranked = p[o] * len(p) / (np.arange(len(p)) + 1)
    q[o] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    return q


def call(gl, floor=FLOOR, alpha=ALPHA, use_floor=True):
    """gl: DataFrame with columns gene, y (guide-level values); NTC rows give the null."""
    ntc = gl.loc[gl.gene == "NTC", "y"].values
    rows = []
    for g, sub in gl[gl.gene != "NTC"].groupby("gene"):
        p = mannwhitneyu(sub.y.values, ntc, alternative="less").pvalue
        rows.append((g, float(np.median(sub.y)), p))
    df = pd.DataFrame(rows, columns=["gene", "lfc", "p"]).set_index("gene")
    df["fdr"] = bh(df.p.values)
    sel = (df.fdr < alpha) & ((df.lfc <= floor) if use_floor else True)
    return df, sorted(df.index[sel])


def frame(y):
    return pd.DataFrame({"gene": lib.gene.values, "y": y.reindex(lib.index).values}, index=lib.index)


LFC = guide_lfc()
CN = gene_cn()
GL = frame(LFC)
GL["log2cn"] = np.log2(GL.gene.map(CN) / 2.0)


# ---- correct variants ------------------------------------------------------
def correct_regression():
    """Robust regression of guide LFC on log2 copy number (gain and loss terms), calls on residuals."""
    t = GL[GL.gene != "NTC"]
    xp, xm, y = t.log2cn.clip(lower=0).values, (-t.log2cn).clip(lower=0).values, t.y.values
    fit = least_squares(lambda b: b[0] + b[1] * xp + b[2] * xm - y, x0=[0, 0, 0], loss="soft_l1", f_scale=0.5)
    a, bp, bm = fit.x
    r = GL.copy()
    xp_all, xm_all = r.log2cn.fillna(0).clip(lower=0), (-r.log2cn.fillna(0)).clip(lower=0)
    r["y"] = r.y - (bp * xp_all + bm * xm_all)
    return call(r)[1], dict(slope_gain=bp, shift_loss=bm)


def correct_segment_median():
    """CRISPRcleanR-like: subtract, per copy-number segment with >=5 genes and |log2(CN/2)|>0.3, the segment's
    median guide LFC (robust to the few true essentials inside a segment)."""
    r = GL.copy()
    gcn = CN.reset_index().rename(columns={"index": "gene"})
    sid = {}
    for _, row in coords.iterrows():
        mid = (row.start + row.end) / 2
        k = segs.index[(segs.chromosome == row.chromosome) & (segs.start <= mid) & (segs.end >= mid)][0]
        sid[row.gene] = k
    r["seg"] = r.gene.map(sid)
    for k, sub in r[r.gene != "NTC"].groupby("seg"):
        cn = segs.loc[k, "copy_number"]
        if sub.gene.nunique() >= 5 and abs(np.log2(cn / 2)) > 0.3:
            r.loc[sub.index, "y"] = sub.y - sub.y.median()
    return call(r)[1]


def correct_exclude():
    """Drop genes in gained regions (CN >= 2.5), then call the rest."""
    keep = set(CN.index[CN < 2.5]) | {"NTC"}
    return call(GL[GL.gene.isin(keep)])[1]


# ---- traps -----------------------------------------------------------------
def trap_naive_lfc():
    """Gene-level LFC <= -1 and FDR<0.05, no copy-number handling."""
    return call(GL)[1]


def trap_rra_like_fdr():
    """MAGeCK-style: every gene significantly depleted at FDR<0.05, no effect-size floor."""
    return call(GL, use_floor=False)[1]


def trap_topN():
    """Top N genes by gene-level LFC, N = number significant at FDR<0.05."""
    df, sig = call(GL, use_floor=False)
    return sorted(df.sort_values("lfc").index[:len(sig)])


def info_cn_fdr_no_floor():
    """CN-corrected but ignoring the stated effect-size floor (documented, not required to pass)."""
    t = GL[GL.gene != "NTC"]
    xp = t.log2cn.clip(lower=0).values
    r = GL.copy()
    fit = least_squares(lambda b: b[0] + b[1] * xp - t.y.values, x0=[0, 0], loss="soft_l1", f_scale=0.5)
    r["y"] = r.y - fit.x[1] * r.log2cn.fillna(0).clip(lower=0)
    return call(r, use_floor=False)[1]


def blindspot_naive_floor_minus2():
    """BLIND SPOT: naive caller with a much stricter floor (-2) than the task states; passes (0.86/0.80)."""
    return call(GL, floor=-2.0)[1]


def blindspot_naive_top40():
    """BLIND SPOT: naive LFC ranking truncated at an arbitrary N=40 (needs luck/knowledge of N); passes."""
    df, _ = call(GL, use_floor=False)
    return sorted(df.sort_values("lfc").index[:40])


def ans(genes):
    return {"essential_genes": genes}


if __name__ == "__main__":
    truth = json.loads((HERE / "truth.json").read_text())
    T = set(truth["essential_genes"])
    amp = set(truth["amplicon_non_essential"])
    weak = set(truth["weak_genes"])
    reg_genes, reg_info = correct_regression()
    print("regression fit:", {k: round(v, 3) for k, v in reg_info.items()},
          "(truth gain slope -1.2*mean(e)~-0.9)")
    ev = HERE / "eval.json"
    solvers = [("correct: regression", reg_genes, True),
               ("correct: segment-median", correct_segment_median(), True),
               ("correct: exclude gained", correct_exclude(), True),
               ("trap: naive LFC<=-1", trap_naive_lfc(), False),
               ("trap: RRA-like FDR only", trap_rra_like_fdr(), False),
               ("trap: top-N by LFC", trap_topN(), False),
               ("info: CN-corr, no floor", info_cn_fdr_no_floor(), None),
               ("BLIND SPOT naive floor -2", blindspot_naive_floor_minus2(), True),
               ("BLIND SPOT naive top-40", blindspot_naive_top40(), True)]
    ok = True
    print(f"{'solver':26s} {'n':>4s} {'TP':>3s} {'FP-amp':>6s} {'FP-weak':>7s} {'FP-oth':>6s} {'prec':>5s} {'rec':>5s} pass")
    for name, g, expect in solvers:
        s = set(g)
        tp = len(s & T)
        r = grade(ev, ans(g))
        print(f"{name:26s} {len(s):4d} {tp:3d} {len(s & amp):6d} {len(s & weak):7d} {len(s - T - amp - weak):6d} "
              f"{tp / max(1, len(s)):5.2f} {tp / len(T):5.2f} {r['passed']}")
        if expect is not None:
            ok &= (r["passed"] == expect)
    print("DISCRIMINATES (as designed)" if ok else "!! DOES NOT DISCRIMINATE")
    sys.exit(0 if ok else 1)
