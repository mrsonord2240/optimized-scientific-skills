"""Genome-equivalent and Poisson *sampling* calculator.

These functions describe whether mutant templates are physically present in an
aliquot. They do not model extraction/recovery, consensus depth, background
error, false positives, or a calling rule and therefore do not estimate an
achieved assay LoD.
"""
# Reference: numpy 1.26+, scipy 1.12+ | Verify API if version differs

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


if __name__ == "__main__":
    print("THEORETICAL POISSON SAMPLING ONLY - NOT AN ACHIEVED ASSAY LoD")
    print("Assumptions: intact input mass, perfect recovery/detection of present templates, no background errors or false positives.")
    print(f"Convention: {DEFAULT_GE_PER_NG:.2f} haploid GE/ng (3.3 pg per haploid genome).")
    print("\nSampling probability for various input masses at 0.1% VAF:")
    for ng in (3.0, 10.0, 30.0):
        ge = genome_equivalents(ng)
        prob = sampling_probability(ng, 0.001)
        print(f"  {ng:>6.1f} ng = {ge:>8.0f} GE | P(template present) = {prob:.3f}")

    print("\nGE for 95% template-presence probability:")
    for vaf in (1e-3, 1e-4, 1e-5):
        ge = ge_for_sampling_probability(vaf)
        print(f"  {vaf:.0e} VAF: {ge:.0f} GE ({ge / DEFAULT_GE_PER_NG:.1f} ng)")
