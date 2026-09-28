"""Panel-integrated LoD calculations.

Combine independent per-locus detection probabilities into the panel-level
detection probability that a bespoke MRD assay actually achieves, and find
the integrated LoD.

Approach: Treat each tracked locus as an independent Poisson sampler at the
same tumor VAF; a panel positive call requires at least k loci detected, so
the panel detection probability is the binomial-tail over the per-locus
probabilities -- this is why summing 16-50 loci reaches ppm.
"""
# Reference: numpy 1.26+, scipy 1.12+ | Verify API if version differs

import numpy as np
from scipy.stats import poisson, binom

GE_PER_NG = 330


def detection_probability(input_ng, vaf, min_mutant_molecules=1):
    """P(at least min_mutant_molecules present) under Poisson(lambda = GE * VAF)."""
    lam = input_ng * GE_PER_NG * vaf
    return float(poisson.sf(min_mutant_molecules - 1, lam))


def panel_detection_probability(input_ng, vaf, n_loci, min_loci_positive=2):
    """P(>= min_loci_positive of n_loci detected); >=2-of-N is the Signatera-style positivity rule.

    Each locus is an independent Poisson sampler at the same tumor VAF.
    """
    per_locus = detection_probability(input_ng, vaf)
    return float(binom.sf(min_loci_positive - 1, n_loci, per_locus))


def panel_integrated_lod95(input_ng, n_loci, min_loci_positive=2, grid=None):
    """Lowest VAF on a log grid where the >=k-of-N panel call hits 95%.

    Searches grid from 1e-6 to 1e-2 by default.
    """
    grid = np.logspace(-6, -2, 400) if grid is None else np.asarray(grid)
    probs = [panel_detection_probability(input_ng, v, n_loci, min_loci_positive) for v in grid]
    hits = grid[np.asarray(probs) >= 0.95]
    return float(hits.min()) if hits.size else float('nan')


if __name__ == '__main__':
    print('Per-locus vs panel-integrated detection at 30 ng, 0.01% VAF:')
    print(f'  Single locus      P(detect) = {detection_probability(30, 1e-4):.3f}')
    for n in (16, 32, 48):
        panel_prob = panel_detection_probability(30, 1e-4, n)
        print(f'  {n:>2}-variant panel  P(>=2 detected) = {panel_prob:.3f}')

    print('\nPanel-integrated LoD95 (VAF at 95% detection):')
    for n in (16, 32, 48):
        lod95 = panel_integrated_lod95(30, n)
        print(f'  {n:>2}-variant panel: {lod95:.2e} VAF')
