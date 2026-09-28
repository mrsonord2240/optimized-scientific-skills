"""Genome equivalents and Poisson detection calculator.

Convert an input mass and target VAF into an expected mutant-molecule count
and a detection probability, so a sensitivity claim is anchored to molecules
rather than to a VAF alone.

Approach: Convert ng to haploid genome equivalents (~330/ng), set
lambda = input_GE x VAF, and read the detection probability as a Poisson
tail P(X >= k) = 1 - cdf(k-1, lambda); invert for the minimum GE that puts
lambda at the >=3 sampling-detection threshold.
"""
# Reference: numpy 1.26+, scipy 1.12+ | Verify API if version differs

import numpy as np
from scipy.stats import poisson

GE_PER_NG = 330  # haploid ~3.3 pg -> strict 1 ng / 3.3 pg = 303; 330 is the common diploid-6.6 pg/rounding convention


def genome_equivalents(input_ng):
    """Convert ng to haploid genome equivalents."""
    return input_ng * GE_PER_NG


def detection_probability(input_ng, vaf, min_mutant_molecules=1):
    """P(at least min_mutant_molecules present) under Poisson(lambda = GE * VAF)."""
    lam = genome_equivalents(input_ng) * vaf
    return float(poisson.sf(min_mutant_molecules - 1, lam))


def ge_for_sampling_detection(vaf, target_lambda=3.0):
    """GE needed so lambda >= 3 -> ~95% chance the mutant molecule is present at all.

    1 - e^-3 = 0.95; ~30,000 GE (~91 ng at 330 GE/ng) for a single 1e-4 variant.
    """
    return target_lambda / vaf


if __name__ == '__main__':
    # A 0.1% variant on 3.0 ng (~990 GE) has lambda ~= 1 -> ~63% detected, ~37% missed by sampling alone.
    print('Sampling ceiling for various input masses at 0.1% VAF:')
    for ng in (3.0, 10.0, 30.0):
        ge = genome_equivalents(ng)
        prob = detection_probability(ng, 0.001)
        print(f'  {ng:>6.1f} ng = {ge:>8.0f} GE | P(detect) = {prob:.3f}')

    print('\nGE for 95% sampling-detection:')
    print(f'  0.1% variant (1e-3): {ge_for_sampling_detection(1e-3):.0f} GE ({ge_for_sampling_detection(1e-3)/GE_PER_NG:.1f} ng)')
    print(f'  0.01% variant (1e-4): {ge_for_sampling_detection(1e-4):.0f} GE ({ge_for_sampling_detection(1e-4)/GE_PER_NG:.1f} ng)')
    print(f'  0.001% variant (1e-5): {ge_for_sampling_detection(1e-5):.0f} GE ({ge_for_sampling_detection(1e-5)/GE_PER_NG:.1f} ng)')
