"""Worked, non-clinical examples for analytical detection-limit calculations."""
# Reference: numpy 1.26+, scipy 1.12+, statsmodels 0.14+ | Verify API if version differs

from numbers import Integral, Real
from pathlib import Path
import sys

import numpy as np
from scipy.stats import norm

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from ge_and_poisson import (  # noqa: E402
    DEFAULT_GE_PER_NG,
    ge_for_sampling_probability,
    genome_equivalents,
    sampling_probability,
)
from lod95_probit import limit_of_blank, lod95_probit  # noqa: E402
from panel_integrated_lod import (  # noqa: E402
    SAMPLING_ONLY_ASSUMPTIONS,
    sampling_only_panel_probability,
    theoretical_sampling_vaf95_lower_bound,
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


def simulate_dilution_series(true_lod_vaf, levels, replicates, seed=0, probit_slope=2.0):
    """Simulate a contrived binary detection curve calibrated at ``true_lod_vaf``.

    The documented teaching model is

    ``P(detect) = Phi(Phi^-1(0.95) + probit_slope * log10(level / true_lod_vaf))``.

    Thus ``true_lod_vaf`` materially sets the 95% point. This is not a Poisson
    template-presence model and is not evidence of assay performance.
    """
    true_lod_vaf = _finite_real("true_lod_vaf", true_lod_vaf, minimum=0.0, maximum=1.0, strict_minimum=True)
    probit_slope = _finite_real("probit_slope", probit_slope, minimum=0.0, strict_minimum=True)
    if isinstance(replicates, (bool, np.bool_)) or not isinstance(replicates, Integral) or replicates < 3:
        raise ValueError("replicates must be an integer of at least 3.")
    try:
        levels = np.asarray(levels, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("levels must contain numeric VAF values.") from exc
    if levels.ndim != 1 or levels.size < 3:
        raise ValueError("levels must be a one-dimensional array with at least 3 VAF values.")
    if not np.all(np.isfinite(levels)) or np.any(levels <= 0.0) or np.any(levels > 1.0):
        raise ValueError("levels must contain only finite VAF values in the interval (0, 1].")
    if np.any(np.diff(levels) <= 0.0):
        raise ValueError("levels must be strictly increasing.")

    probabilities = norm.cdf(norm.ppf(0.95) + probit_slope * np.log10(levels / true_lod_vaf))
    rng = np.random.default_rng(seed)
    vaf_out = np.repeat(levels, replicates)
    probability_out = np.repeat(probabilities, replicates)
    detected = rng.binomial(1, probability_out)
    return vaf_out, detected


def main():
    print("THEORETICAL POISSON SAMPLING ONLY - NOT AN ACHIEVED ASSAY LoD")
    print(f"Convention: {DEFAULT_GE_PER_NG:.2f} haploid GE/ng (3.3 pg per haploid genome).")
    for ng in (3.3, 10, 30):
        print(
            f"  {ng:>5} ng = {genome_equivalents(ng):>6.0f} GE | "
            f"P(template present at 0.1% VAF) = {sampling_probability(ng, 0.001):.3f}"
        )
    required_ge = ge_for_sampling_probability(1e-4)
    print(f"  GE for 95% template presence at 0.01% VAF: {required_ge:.0f} ({required_ge / DEFAULT_GE_PER_NG:.0f} ng)")
    print("  Excludes recovery, consensus depth, background error, false positives, and calling rules.")

    print("\nTHEORETICAL MULTI-LOCUS SAMPLING LOWER BOUND - NOT PANEL LoD95")
    for assumption in SAMPLING_ONLY_ASSUMPTIONS:
        print(f"  assumption: {assumption}")
    for n_loci in (16, 48):
        probability = sampling_only_panel_probability(30, 1e-4, n_loci)
        result = theoretical_sampling_vaf95_lower_bound(30, n_loci)
        print(
            f"  {n_loci:>2} loci: P(>=2 templates present at 1e-4)={probability:.3f}; "
            f"sampling-only VAF lower bound={result['vaf_lower_bound']:.2e}"
        )

    print("\nCONTRIVED PROBIT TEACHING MODEL - NOT VALIDATION EVIDENCE")
    levels = [2e-4, 5e-4, 1e-3, 2e-3, 5e-3]
    vaf, detected = simulate_dilution_series(true_lod_vaf=2e-3, levels=levels, replicates=100, seed=7)
    blanks = np.clip(np.random.default_rng(1).normal(0.0002, 0.00005, size=20), 0.0, None)
    result = lod95_probit(vaf, detected)
    lower, upper = result["lod95_confidence_interval"]
    print(f"  LoB (Gaussian shortcut) = {limit_of_blank(blanks):.6f} VAF")
    print(f"  fitted LoD95 = {result['lod95_vaf']:.6f} VAF; 95% CI [{lower:.6f}, {upper:.6f}]")
    print(
        f"  diagnostics: converged={result['converged']}, slope={result['slope']:.3f}, "
        f"n={result['n_observations']}, brackets_0.95={result['brackets_0_95']}"
    )


if __name__ == "__main__":
    main()
