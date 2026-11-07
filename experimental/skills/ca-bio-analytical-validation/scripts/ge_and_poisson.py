"""Genome-equivalent and Poisson *sampling* calculator.

These functions describe whether mutant templates are physically present in an
aliquot. They do not model extraction/recovery, consensus depth, background
error, false positives, or a calling rule and therefore do not estimate an
achieved assay LoD.
"""

from numbers import Integral, Real

import numpy as np
from scipy.optimize import brentq
from scipy.stats import poisson

# 1 ng / 3.3 pg per haploid genome = 303.03 haploid genome equivalents.
DEFAULT_GE_PER_NG = 1000.0 / 3.3


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


def genome_equivalents(input_ng, ge_per_ng=DEFAULT_GE_PER_NG):
    """Convert positive input mass to haploid genome equivalents.

    ``ge_per_ng`` is explicit because laboratories may adopt a documented
    convention other than the 3.3-pg default.
    """
    input_ng = _finite_real("input_ng", input_ng, minimum=0.0, strict_minimum=True)
    ge_per_ng = _finite_real("ge_per_ng", ge_per_ng, minimum=0.0, strict_minimum=True)
    return input_ng * ge_per_ng


def sampling_probability(input_ng, vaf, min_mutant_molecules=1, ge_per_ng=DEFAULT_GE_PER_NG):
    """Return P(>=k mutant templates present), not assay sensitivity."""
    vaf = _finite_real("vaf", vaf, minimum=0.0, maximum=1.0, strict_minimum=True)
    min_mutant_molecules = _positive_integer("min_mutant_molecules", min_mutant_molecules)
    lam = genome_equivalents(input_ng, ge_per_ng=ge_per_ng) * vaf
    return float(poisson.sf(min_mutant_molecules - 1, lam))


def ge_for_sampling_probability(vaf, target_probability=0.95, min_mutant_molecules=1):
    """Return GE required for a target template-presence probability.

    This is an ideal sampling requirement only. It excludes recovery, read
    support, background error, false positives, and assay calling behavior.
    """
    vaf = _finite_real("vaf", vaf, minimum=0.0, maximum=1.0, strict_minimum=True)
    target_probability = _finite_real(
        "target_probability", target_probability, minimum=0.0, maximum=1.0, strict_minimum=True
    )
    if target_probability >= 1.0:
        raise ValueError("target_probability must be less than 1.0.")
    min_mutant_molecules = _positive_integer("min_mutant_molecules", min_mutant_molecules)

    upper = float(max(1, min_mutant_molecules))
    while poisson.sf(min_mutant_molecules - 1, upper) < target_probability:
        upper *= 2.0
    required_lambda = float(
        brentq(
            lambda lam: poisson.sf(min_mutant_molecules - 1, lam) - target_probability,
            np.finfo(float).tiny,
            upper,
        )
    )
    return required_lambda / vaf


def main():
    import argparse
    import csv
    import sys

    ap = argparse.ArgumentParser(description="Genome equivalents and Poisson template-presence probability (sampling only).")
    ap.add_argument("--input-ng", type=float, nargs="+", required=True, help="cfDNA input mass in ng (one or more)")
    ap.add_argument("--vaf", type=float, nargs="+", required=True, help="variant allele fraction as a fraction, 0.1%% = 0.001 (one or more)")
    ap.add_argument("--min-molecules", type=int, default=1, help="mutant templates required to count as present")
    ap.add_argument("--target", type=float, default=0.95, help="target probability for the required-GE column")
    ap.add_argument("--ge-per-ng", type=float, default=DEFAULT_GE_PER_NG, help="state the lab convention; default 1000/3.3")
    ap.add_argument("--out", default="ge_poisson_results.tsv")
    a = ap.parse_args()
    if any(1 < v <= 100 for v in a.vaf):
        sys.exit("REFUSED: --vaf is a fraction, not a percent (0.1% is 0.001). Values above 1 look like percents.")
    rows = []
    try:
        for ng in a.input_ng:
            for vaf in a.vaf:
                ge = genome_equivalents(ng, ge_per_ng=a.ge_per_ng)
                need = ge_for_sampling_probability(vaf, a.target, a.min_molecules)
                rows.append({
                    "input_ng": ng, "ge_per_ng": a.ge_per_ng, "genome_equivalents": round(ge, 2), "vaf": vaf,
                    "expected_mutant_molecules": round(ge * vaf, 4), "min_molecules": a.min_molecules,
                    "p_template_present": round(sampling_probability(ng, vaf, a.min_molecules, a.ge_per_ng), 6),
                    "target_probability": a.target, "ge_required": round(need, 1),
                    "ng_required": round(need / a.ge_per_ng, 2),
                    "scope": "sampling only; excludes recovery, consensus depth, background error, calling rule"})
    except ValueError as exc:
        sys.exit(f"REFUSED: {exc}")
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
