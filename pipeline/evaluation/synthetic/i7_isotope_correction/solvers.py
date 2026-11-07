"""Reference solvers for synthetic task i7_isotope_correction.

`correct()`  : full MID correction. Response matrix from the MEASURED ION formula (natural
               C on every non-labelled carbon incl. derivatisation carbons, plus H, N, O, Si),
               tracer purity applied at labelled positions, non-negative least squares over all
               reported isotopologues, enrichment = sum i*MID_i / backbone_carbons.
`trap_*()`   : shortcuts a person might take. Each must FAIL the graders.

Usage: python solvers.py          (prints table; exit 0 iff DISCRIMINATES)
"""
import json
import re
import sys
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from grading import grade  # noqa: E402

D = HERE / "data"
AREAS = pd.read_csv(D / "isotopologue_areas.tsv", sep="\t")
MET = pd.read_csv(D / "metabolites.tsv", sep="\t").set_index("metabolite")
MCOLS = [c for c in AREAS.columns if re.fullmatch(r"M\d+", c)]
EV = json.loads((HERE / "eval.json").read_text(encoding="utf-8"))
TRUTH = EV["graders"][0]["config"]["ground_truth"]
REQ = [(k.rsplit("_", 1)[0], k.rsplit("_", 1)[1]) for k in TRUTH]          # (metabolite, "4h")

# IUPAC natural abundances (mass shift 0,1,2)
NAT = {"C": [0.9893, 0.0107], "H": [0.999885, 0.000115], "N": [0.99636, 0.00364],
       "O": [0.99757, 0.00038, 0.00205], "Si": [0.92223, 0.04685, 0.03092]}


def counts(formula):
    return {e: int(n or 1) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}


def poly_pow(p, k, L):
    out = np.zeros(L)
    out[0] = 1
    for _ in range(k):
        out = np.convolve(out, p)[:L]
    return out


def response(formula, n, purity, L, elements=("C", "H", "N", "O", "Si"), n_carbon_label=None):
    """(L, n+1) matrix: column k = measured pattern of a molecule with k tracer-derived carbons."""
    el = counts(formula)
    A = np.zeros((L, n + 1))
    for k in range(n + 1):
        lab = np.zeros(L)
        for j in range(min(k, L - 1) + 1):
            lab[j] = comb(k, j) * purity ** j * (1 - purity) ** (k - j)
        v = np.convolve(lab, poly_pow(np.array(NAT["C"]), el["C"] - k, L))[:L] if "C" in elements else lab
        for e in ("H", "N", "O", "Si"):
            if e in elements and e in el:
                v = np.convolve(v, poly_pow(np.array(NAT[e]), el[e], L))[:L]
        A[:, k] = v
    return A


def raw_vec(metab, sample):
    r = AREAS[(AREAS.metabolite == metab) & (AREAS["sample"] == sample)].iloc[0]
    v = r[MCOLS].astype(float).values
    return v[~np.isnan(v)]


def correct_mid(metab, sample, purity=None, formula=None, elements=("C", "H", "N", "O", "Si")):
    m = MET.loc[metab]
    n = int(m.backbone_carbons)
    purity = float(m.tracer_13C_purity) if purity is None else purity
    formula = m.ion_formula if formula is None else formula
    y = raw_vec(metab, sample)
    A = response(formula, n, purity, len(y), elements)
    c, _ = nnls(A, y)
    return c / c.sum()


def enrich(mid, n):
    return float((np.arange(len(mid)) * mid).sum() / n)


def answer_from(fn):
    return {f"{m}_{t}": round(fn(m, f"Glc13C_{t}"), 4) for m, t in REQ}


# ---------------------------------------------------------------- solvers
def correct():
    return answer_from(lambda m, s: enrich(correct_mid(m, s), int(MET.loc[m].backbone_carbons)))


def trap_raw():
    """Uncorrected: normalise reported peak areas, sum i*f_i / n."""
    def f(m, s):
        v = raw_vec(m, s)
        return enrich(v / v.sum(), int(MET.loc[m].backbone_carbons))
    return answer_from(f)


def trap_carbon_only():
    """Correct natural 13C only (all carbons of the measured ion), ignoring H/N/O/Si."""
    return answer_from(lambda m, s: enrich(correct_mid(m, s, elements=("C",)), int(MET.loc[m].backbone_carbons)))


def trap_backbone_formula():
    """Correct with the underivatised metabolite formula, ignoring the derivatisation atoms in the ion."""
    UNDERIV = {"lactate": "C3H5O3", "alanine": "C3H6NO2", "citrate": "C6H7O7", "glutamate": "C5H8NO4"}
    return answer_from(lambda m, s: enrich(correct_mid(m, s, formula=UNDERIV[m]), int(MET.loc[m].backbone_carbons)))


def trap_no_purity():
    """Full natural-abundance correction, tracer assumed 100% pure."""
    return answer_from(lambda m, s: enrich(correct_mid(m, s, purity=1.0), int(MET.loc[m].backbone_carbons)))


def trap_subtract_control_enrichment():
    """Raw enrichment minus raw enrichment of the unlabelled control."""
    def f(m, s):
        n = int(MET.loc[m].backbone_carbons)
        e = lambda smp: enrich(raw_vec(m, smp) / raw_vec(m, smp).sum(), n)
        return e(s) - e("Glc12C_24h")
    return answer_from(f)


def trap_subtract_control_mid():
    """Subtract control's M1..Mn fractions from the sample's, put the remainder in M0."""
    def f(m, s):
        n = int(MET.loc[m].backbone_carbons)
        a, b = raw_vec(m, s), raw_vec(m, "Glc12C_24h")
        sub = np.clip(a[1:n + 1] / a.sum() - b[1:n + 1] / b.sum(), 0, None)
        return enrich(np.concatenate([[1 - sub.sum()], sub]), n)
    return answer_from(f)


def trap_wrong_n():
    """Correct MID but divide by carbons in the measured ion formula instead of backbone carbons."""
    def f(m, s):
        c = correct_mid(m, s)
        return enrich(c, counts(MET.loc[m].ion_formula)["C"])
    return answer_from(f)


def isocor_check():
    """Independent cross-check with the isocor package (pip install isocor). Not gating."""
    from isocor.mscorrectors import LowResMetaboliteCorrector as C
    cfg = {"lactate": ("C3H4O3", "C8H21Si2"), "alanine": ("C3H5NO2", "C8H21Si2"),
           "citrate": ("C6H7O7", ""), "glutamate": ("C5H6NO4", "C14H36Si3")}   # moiety, derivative
    p = float(MET.tracer_13C_purity.iloc[0])

    def f(m, s):
        n = int(MET.loc[m].backbone_carbons)
        c = C(cfg[m][0], "13C", data_isotopes=C.DEFAULT_ISODATA, tracer_purity=[1 - p, p],
              correct_NA_tracer=True, derivative_formula=cfg[m][1])
        return float(c.correct(raw_vec(m, s)[: n + 1])[3])
    return answer_from(f)


SOLVERS = [("correct (should PASS)", correct, True),
           ("trap raw MID", trap_raw, False),
           ("trap 13C-only correction", trap_carbon_only, False),
           ("trap backbone formula only", trap_backbone_formula, False),
           ("trap ignore tracer purity", trap_no_purity, False),
           ("trap subtract ctrl enrichment", trap_subtract_control_enrichment, False),
           ("trap subtract ctrl MID", trap_subtract_control_mid, False),
           ("trap n = all ion carbons", trap_wrong_n, False)]


def main():
    ev = HERE / "eval.json"
    ok = True
    print("truth:", TRUTH)
    for name, fn, expect in SOLVERS:
        a = fn()
        r = grade(ev, a)
        err = {k: a[k] - TRUTH[k] for k in TRUTH}
        nfail = sum(abs(v) > 0.01 for v in err.values())
        print(f"{name:30s} all={r['passed']!s:5s} fields_failing={nfail}/{len(err)} max|err|={max(abs(v) for v in err.values()):.4f} "
              f"min|err|={min(abs(v) for v in err.values()):.4f}")
        print("     ", {k: f"{v:+.3f}" for k, v in err.items()})
        label = "" if r["passed"] == expect else "  <-- UNEXPECTED"
        if r["passed"] and not expect:
            label = "  <-- BLIND SPOT"
        if label:
            print("     ", label)
        ok &= (r["passed"] == expect)
    try:
        a = isocor_check()
        print(f"isocor cross-check (not gating) all={grade(ev, a)['passed']} max|err|={max(abs(a[k] - TRUTH[k]) for k in TRUTH):.4f}")
    except Exception as e:  # isocor not installed
        print("isocor cross-check skipped:", type(e).__name__)
    print("DISCRIMINATES (as designed)" if ok else "!! DOES NOT DISCRIMINATE")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
