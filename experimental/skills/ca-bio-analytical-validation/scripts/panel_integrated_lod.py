"""Idealized multi-locus Poisson sampling calculations.

The outputs are theoretical lower bounds on an assay's attainable VAF
threshold, never achieved panel LoD95 values. The model assumes equal-VAF,
independent loci and perfect recovery/calling with no background errors.
"""

from numbers import Integral, Real

import numpy as np
from scipy.stats import binom, poisson

DEFAULT_GE_PER_NG = 1000.0 / 3.3
SAMPLING_ONLY_ASSUMPTIONS = (
    "equal VAF at every tracked locus",
    "independent locus sampling",
    "perfect molecular recovery and detection of every present mutant template",
    "no consensus-depth failures, background errors, false positives, or empirical calling effects",
)


def _finite_real(name, value, *, minimum=None, maximum=None, strict_minimum=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not np.isfinite(value):
        raise ValueError(f"{name} must be a finite real number.")
    value = float(value)
    if minimum is not None and (value <= minimum if strict_minimum else value < minimum):
        relation = "greater than" if strict_minimum else "at least"
        raise ValueError(f"{name} must be {relation} {minimum}.")
    if maximum is not None and value > maximum:
        raise ValueError(f"{name} must be at most {maximum}.")
    return value


def _positive_integer(name, value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 1:
        raise ValueError(f"{name} must be an integer of at least 1.")
    return int(value)


def _panel_rule(n_loci, min_loci_positive):
    n_loci = _positive_integer("n_loci", n_loci)
    min_loci_positive = _positive_integer("min_loci_positive", min_loci_positive)
    if min_loci_positive > n_loci:
        raise ValueError("min_loci_positive must not exceed n_loci.")
    return n_loci, min_loci_positive


def sampling_only_locus_probability(input_ng, vaf, ge_per_ng=DEFAULT_GE_PER_NG):
    """P(at least one mutant template is present); not assay sensitivity."""
    input_ng = _finite_real("input_ng", input_ng, minimum=0.0, strict_minimum=True)
    vaf = _finite_real("vaf", vaf, minimum=0.0, maximum=1.0, strict_minimum=True)
    ge_per_ng = _finite_real("ge_per_ng", ge_per_ng, minimum=0.0, strict_minimum=True)
    return float(poisson.sf(0, input_ng * ge_per_ng * vaf))


def sampling_only_panel_probability(input_ng, vaf, n_loci, min_loci_positive=2, ge_per_ng=DEFAULT_GE_PER_NG):
    """Ideal P(>=k of N loci have a mutant template physically present)."""
    n_loci, min_loci_positive = _panel_rule(n_loci, min_loci_positive)
    per_locus = sampling_only_locus_probability(input_ng, vaf, ge_per_ng=ge_per_ng)
    return float(binom.sf(min_loci_positive - 1, n_loci, per_locus))


def theoretical_sampling_vaf95_lower_bound(
    input_ng,
    n_loci,
    min_loci_positive=2,
    grid=None,
    target_probability=0.95,
    ge_per_ng=DEFAULT_GE_PER_NG,
):
    """Return a structured sampling-only VAF lower bound.

    The result is the lowest supplied grid value whose idealized k-of-N
    template-presence probability reaches ``target_probability``. It must not
    be reported as an assay LoD95.
    """
    input_ng = _finite_real("input_ng", input_ng, minimum=0.0, strict_minimum=True)
    n_loci, min_loci_positive = _panel_rule(n_loci, min_loci_positive)
    target_probability = _finite_real(
        "target_probability", target_probability, minimum=0.0, maximum=1.0, strict_minimum=True
    )
    if target_probability >= 1.0:
        raise ValueError("target_probability must be less than 1.0.")
    ge_per_ng = _finite_real("ge_per_ng", ge_per_ng, minimum=0.0, strict_minimum=True)

    values = np.logspace(-6, -2, 400) if grid is None else np.asarray(grid, dtype=float)
    if values.ndim != 1 or values.size < 2:
        raise ValueError("grid must be a one-dimensional array with at least 2 VAF values.")
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0) or np.any(values > 1.0):
        raise ValueError("grid must contain only finite VAF values in the interval (0, 1].")
    if np.any(np.diff(values) <= 0.0):
        raise ValueError("grid must be strictly increasing.")

    probs = np.asarray(
        [sampling_only_panel_probability(input_ng, vaf, n_loci, min_loci_positive, ge_per_ng) for vaf in values]
    )
    hits = np.flatnonzero(probs >= target_probability)
    if hits.size == 0:
        raise ValueError("grid does not reach target_probability; extend its upper VAF bound.")
    index = int(hits[0])
    return {
        "vaf_lower_bound": float(values[index]),
        "sampling_probability": float(probs[index]),
        "target_probability": target_probability,
        "input_ng": input_ng,
        "ge_per_ng": ge_per_ng,
        "n_loci": n_loci,
        "min_loci_positive": min_loci_positive,
        "assumptions": list(SAMPLING_ONLY_ASSUMPTIONS),
        "interpretation": "theoretical sampling-only lower bound; not an achieved assay LoD95",
    }


def main():
    import argparse
    import csv
    import sys

    ap = argparse.ArgumentParser(description="Theoretical multi-locus SAMPLING lower bound; not an assay LoD95.")
    ap.add_argument("--input-ng", type=float, required=True)
    ap.add_argument("--n-loci", type=int, nargs="+", required=True, help="panel size(s)")
    ap.add_argument("--min-loci", type=int, default=2, help="loci that must have a template present (k of N)")
    ap.add_argument("--target", type=float, default=0.95)
    ap.add_argument("--ge-per-ng", type=float, default=DEFAULT_GE_PER_NG)
    ap.add_argument("--out", default="panel_bound_results.tsv")
    a = ap.parse_args()
    rows = []
    try:
        for n in a.n_loci:
            r = theoretical_sampling_vaf95_lower_bound(a.input_ng, n, a.min_loci, target_probability=a.target, ge_per_ng=a.ge_per_ng)
            rows.append({"input_ng": a.input_ng, "ge_per_ng": a.ge_per_ng, "n_loci": n, "min_loci_positive": a.min_loci,
                         "target_probability": a.target, "vaf_sampling_lower_bound": r["vaf_lower_bound"],
                         "sampling_probability": round(r["sampling_probability"], 5),
                         "interpretation": r["interpretation"]})
    except ValueError as exc:
        sys.exit(f"REFUSED: {exc}")
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("THEORETICAL SAMPLING LOWER BOUND, NOT ASSAY LoD95. Assumptions:")
    for text in SAMPLING_ONLY_ASSUMPTIONS:
        print("  -", text)
    for r in rows:
        print(r)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
