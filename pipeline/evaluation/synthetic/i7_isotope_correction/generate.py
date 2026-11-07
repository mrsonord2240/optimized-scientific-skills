"""Synthetic eval task i7_isotope_correction.

Everything here is FABRICATED. Design source: the published pitfall that raw
mass-isotopologue distributions (MIDs) contain natural-abundance isotopes of all
elements in the measured ion (and of derivatisation atoms) plus tracer impurity, so
enrichment claims need a full correction (Buescher et al. 2015 Curr Opin Biotechnol;
Midani et al. 2017).

World
  * [U-13C6]-glucose tracer, 98 atom% 13C at each labelled position.
  * Four metabolites; the measured ion formula is given (some are TBDMS GC-MS
    fragments containing Si and derivatisation carbons that can never be labelled).
  * True corrected MIDs (distribution of the number of tracer-derived carbons) are set
    explicitly for 4 time points; a culture on unlabelled glucose is the control.
  * Measured MID = sum_k c_k * [Binom(k, purity) 13C at labelled positions
                     (x) natural 13C on all other carbons (x) natural H,N,O,Si].
    Peaks are truncated at M+(n+2), scaled by a random per-sample factor, and given
    ~0.5% relative Gaussian noise.

Run: python generate.py -> data/, truth.json, eval.json
"""
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)
rng = np.random.default_rng(20260919)

PURITY = 0.98
TOL = 0.01
NOISE = 0.005
EXTRA = 2  # measured isotopologues beyond M+n

# IUPAC-style natural abundances, by mass shift 0,1,2 relative to the lightest isotope.
NAT = {
    "C": [0.9893, 0.0107],
    "H": [0.999885, 0.000115],
    "N": [0.99636, 0.00364],
    "O": [0.99757, 0.00038, 0.00205],
    "Si": [0.92223, 0.04685, 0.03092],
}


def parse_formula(f):
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", f)}


def conv_pow(pmf, k, L):
    out = np.zeros(L)
    out[0] = 1.0
    for _ in range(k):
        out = np.convolve(out, pmf)[:L]
    return out


def binom_pmf(k, p, L):
    from math import comb
    out = np.zeros(L)
    for j in range(min(k, L - 1) + 1):
        out[j] = comb(k, j) * p ** j * (1 - p) ** (k - j)
    return out


def column(k, formula, L, purity=PURITY):
    """Measured pattern (rows M0..M(L-1)) of a molecule with k tracer-derived carbons."""
    el = parse_formula(formula)
    v = binom_pmf(k, purity, L)
    v = np.convolve(v, conv_pow(np.array(NAT["C"]), el["C"] - k, L))[:L]
    for e in ("H", "N", "O", "Si"):
        if e in el:
            v = np.convolve(v, conv_pow(np.array(NAT[e]), el[e], L))[:L]
    return v


# name, measured ion, ion formula, backbone carbons, true MIDs by time (h)
METS = [
    ("lactate", "lactate 2TBDMS [M-57]+, m/z 261", "C11H25O3Si2", 3, {
        1: [0.905, 0.004, 0.006, 0.085],
        4: [0.800, 0.006, 0.014, 0.180],
        8: [0.450, 0.010, 0.030, 0.510],
        24: [0.220, 0.010, 0.030, 0.740]}),
    ("alanine", "alanine 2TBDMS [M-57]+, m/z 260", "C11H26NO2Si2", 3, {
        1: [0.975, 0.002, 0.003, 0.020],
        4: [0.925, 0.005, 0.010, 0.060],
        8: [0.870, 0.008, 0.017, 0.105],
        24: [0.760, 0.012, 0.028, 0.200]}),
    ("citrate", "citrate [M-H]-, m/z 191", "C6H7O7", 6, {
        1: [0.920, 0.004, 0.060, 0.006, 0.008, 0.001, 0.001],
        4: [0.800, 0.008, 0.130, 0.015, 0.035, 0.005, 0.007],
        8: [0.700, 0.010, 0.190, 0.020, 0.055, 0.010, 0.015],
        24: [0.560, 0.012, 0.260, 0.030, 0.090, 0.020, 0.028]}),
    ("glutamate", "glutamate 3TBDMS [M-57]+, m/z 432", "C19H42NO4Si3", 5, {
        1: [0.950, 0.003, 0.038, 0.003, 0.005, 0.001],
        4: [0.860, 0.006, 0.095, 0.008, 0.025, 0.006],
        8: [0.770, 0.008, 0.150, 0.012, 0.045, 0.015],
        24: [0.640, 0.010, 0.220, 0.020, 0.080, 0.030]}),
]
TIMES = [1, 4, 8, 24]
SAMPLES = [(f"Glc13C_{t}h", t) for t in TIMES] + [("Glc12C_24h", 24)]
MAXM = max(n for _, _, _, n, _ in METS) + EXTRA  # 8

# requested combinations (metabolite, time)
REQUESTS = [("lactate", 4), ("lactate", 24), ("alanine", 8), ("citrate", 8), ("glutamate", 4), ("glutamate", 24)]


def enrichment(mid):
    mid = np.asarray(mid, float)
    return float((np.arange(len(mid)) * mid).sum() / (len(mid) - 1))


def main():
    rows, truth = [], {"purity": PURITY, "mids": {}, "enrichment": {}}
    for name, ion, formula, n, mids in METS:
        L = n + EXTRA + 1
        A = np.array([column(k, formula, L) for k in range(n + 1)]).T  # (L, n+1)
        for sample, t in SAMPLES:
            ctrl = sample.startswith("Glc12C")
            c = np.array([1.0] + [0.0] * n) if ctrl else np.array(mids[t])
            assert abs(c.sum() - 1) < 1e-9 and len(c) == n + 1, (name, t)
            if ctrl:
                # unlabelled culture: no tracer-derived carbon; no purity effect (k=0)
                pass
            meas = A @ c
            scale = float(np.exp(rng.uniform(np.log(2e5), np.log(5e7))))
            areas = scale * meas * (1 + rng.normal(0, NOISE, L))
            rows.append((name, sample, t, [f"{a:.0f}" for a in areas] + ["NA"] * (MAXM + 1 - L)))
            truth["mids"][f"{name}|{sample}"] = c.tolist()
            truth["enrichment"][f"{name}|{sample}"] = enrichment(c)

    with open(DATA / "isotopologue_areas.tsv", "w", newline="\n") as f:
        f.write("metabolite\tsample\ttime_h\t" + "\t".join(f"M{i}" for i in range(MAXM + 1)) + "\n")
        for name, sample, t, vals in rows:
            f.write(f"{name}\t{sample}\t{t}\t" + "\t".join(vals) + "\n")

    with open(DATA / "metabolites.tsv", "w", newline="\n") as f:
        f.write("metabolite\tmeasured_ion\tion_formula\tbackbone_carbons\ttracer\ttracer_13C_purity\n")
        for name, ion, formula, n, _ in METS:
            f.write(f"{name}\t{ion}\t{formula}\t{n}\t[U-13C6]-glucose\t{PURITY}\n")

    gt = {f"{m}_{t}h": round(truth["enrichment"][f"{m}|Glc13C_{t}h"], 4) for m, t in REQUESTS}
    truth["requested"] = gt
    (HERE / "truth.json").write_text(json.dumps(truth, indent=2), encoding="utf-8")

    schema = "{" + ", ".join(f'"{k}": <number>' for k in gt) + "}"
    task = (
        "We ran a [U-13C6]-glucose tracing experiment in cultured cells and harvested at 1, 4, 8 and 24 h "
        "(samples Glc13C_1h ... Glc13C_24h). A parallel culture grown on unlabelled glucose and harvested at 24 h "
        "(Glc12C_24h) is also included. Four metabolites were measured by MS: peak areas for each mass isotopologue "
        "(M0, M1, ... = ion mass shifted by 0, 1, ... x 1.00335 Da; NA = not measured) are in data/isotopologue_areas.tsv. "
        "data/metabolites.tsv gives, for each metabolite, the formula of the ion that was actually measured "
        "(this includes any derivatisation atoms retained in the fragment), backbone_carbons (the number of carbons "
        "that come from the metabolite itself, i.e. the only positions that tracer carbon can occupy), the tracer and "
        "tracer_13C_purity (atom fraction of 13C at each labelled position of the tracer). Peak areas are on arbitrary, "
        "sample-specific scales.\n\n"
        "For each of the following metabolite x time-point combinations, report the fractional 13C enrichment "
        "attributable to the tracer: the mean number of tracer-derived carbon atoms per molecule divided by "
        "backbone_carbons, i.e. sum_i i*MID_i / n, where MID_i (i = 0..n, summing to 1) is the fraction of the metabolite "
        "pool carrying i tracer-derived carbons. Report it as a fraction between 0 and 1 (not a percentage), to 3 decimals.\n\n"
        + "; ".join(f"{m} at {t} h" for m, t in REQUESTS)
        + ".\n\nReturn EXACTLY:\n\n<EVAL_ANSWER>\n" + schema + "\n</EVAL_ANSWER>"
    )

    ev = {
        "id": "i7_isotope_correction",
        "replicates": None,
        "synthetic": True,
        "task": task,
        "graders": [{"type": "numeric_tolerance", "config": {
            "tolerances": {k: {"type": "absolute", "lower": TOL, "upper": TOL} for k in gt},
            "ground_truth": gt}}],
        "data_dir": "data",
        "relevant_skills": [
            "bio-metabolomics-isotope-tracing", "bio-metabolomics-targeted-analysis",
            "bio-workflows-metabolomics-pipeline", "bio-metabolomics-normalization-qc"],
        "design_source": "Raw MIDs include natural-abundance isotopes of every element in the measured ion (incl. derivative atoms) and tracer impurity; enrichment needs full correction (Buescher et al. 2015 Curr Opin Biotechnol; Midani et al. 2017).",
    }
    (HERE / "eval.json").write_text(json.dumps(ev, indent=2), encoding="utf-8")
    print("wrote data/, truth.json, eval.json; truth:", gt)


if __name__ == "__main__":
    main()
