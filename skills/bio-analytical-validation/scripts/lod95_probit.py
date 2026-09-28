"""Limit of blank and LoD95 probit fit.

Estimate the VAF at which the assay detects 95% of the time, from a contrived
dilution series, and anchor it to the blank-derived false-positive floor.

Approach: Compute LoB from blank replicates (mean + 1.645*SD, one-sided 95th
pct), then fit a probit GLM of binary detection on log10(VAF) (CLSI EP17 fits
on log concentration) across the dilution series and invert it for the 95%
detection point. The series must bracket the 0.95 crossing -- all-detected
upper levels cause near-complete separation and an unstable slope.
"""
# Reference: numpy 1.26+, scipy 1.12+, statsmodels 0.14+ | Verify API if version differs

import numpy as np
import statsmodels.api as sm
from scipy.stats import norm


def limit_of_blank(blank_signals):
    """LoB = mean + 1.645*SD; one-sided 95th percentile of analyte-free blanks.

    Per CLSI EP17-A2.
    """
    blank_signals = np.asarray(blank_signals, dtype=float)
    return blank_signals.mean() + 1.645 * blank_signals.std(ddof=1)


def lod95_probit(vaf_levels, detected):
    """Probit fit of detection (0/1) on log10(VAF).

    Returns the VAF where P(detect) = 0.95.

    The dilution series must bracket the 0.95 crossing: all-detected upper
    levels cause near-complete separation and an unstable slope.

    CLSI EP17-A2 fits on log concentration to respect the saturating detection
    curve.
    """
    log_vaf = np.log10(np.asarray(vaf_levels, dtype=float))
    y = np.asarray(detected, dtype=float)
    X = sm.add_constant(log_vaf)
    fit = sm.GLM(y, X, family=sm.families.Binomial(link=sm.families.links.Probit())).fit()
    intercept, slope = fit.params
    return 10 ** ((norm.ppf(0.95) - intercept) / slope)


if __name__ == '__main__':
    # Simulated example
    np.random.seed(42)

    # Blank replicates (background error signal)
    blanks = np.random.normal(0.0002, 0.00005, size=20)
    print(f'LoB (blank 95th pct) = {limit_of_blank(blanks):.6f} VAF')

    # Dilution series: VAF levels and detection outcomes
    vaf_levels = [5e-4, 1e-3, 2e-3, 5e-3, 1e-2]
    detected = [0, 1, 1, 1, 1]  # simplified

    try:
        lod95 = lod95_probit(vaf_levels, detected)
        print(f'LoD95 (probit) = {lod95:.6f} VAF')
    except Exception as e:
        print(f'Failed to fit probit: {e}')
