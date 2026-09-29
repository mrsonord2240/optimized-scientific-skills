"""Validated limit-of-blank and replicated LoD95 probit calculations."""
# Reference: numpy 1.26+, scipy 1.12+, statsmodels 0.14+ | Verify API if version differs

import warnings
from numbers import Real

import numpy as np
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.tools.sm_exceptions import PerfectSeparationError, PerfectSeparationWarning


def _one_dimensional_finite(name, values, *, minimum=None, maximum=None):
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain numeric values.") from exc
    if array.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values.")
    if minimum is not None and np.any(array < minimum):
        raise ValueError(f"{name} values must be at least {minimum}.")
    if maximum is not None and np.any(array > maximum):
        raise ValueError(f"{name} values must be at most {maximum}.")
    return array


def limit_of_blank(blank_signals):
    """Gaussian LoB shortcut: mean + 1.645 sample SD.

    At least two finite VAF-scale blank replicates in [0, 1] are required.
    A full EP17 study may require finite-sample or nonparametric treatment.
    """
    blank_signals = _one_dimensional_finite("blank_signals", blank_signals, minimum=0.0, maximum=1.0)
    if blank_signals.size < 2:
        raise ValueError("blank_signals must contain at least 2 replicates.")
    return float(blank_signals.mean() + 1.645 * blank_signals.std(ddof=1))


def lod95_probit(vaf_levels, detected, *, confidence_level=0.95, min_replicates_per_level=3):
    """Fit detection on log10(VAF) and return LoD95, CI, and diagnostics.

    The input must include at least three replicated VAF levels and empirical
    detection rates on both sides of 0.95. Separated, non-converged,
    non-increasing, or otherwise unidentified fits raise stable ``ValueError``
    messages instead of returning a precise-looking result.
    """
    vaf = _one_dimensional_finite("vaf_levels", vaf_levels, minimum=0.0, maximum=1.0)
    y = _one_dimensional_finite("detected", detected, minimum=0.0, maximum=1.0)
    if vaf.size != y.size:
        raise ValueError("vaf_levels and detected must have the same length.")
    if vaf.size == 0:
        raise ValueError("vaf_levels and detected must not be empty.")
    if np.any(vaf <= 0.0):
        raise ValueError("vaf_levels values must be greater than 0 and at most 1.")
    if not np.all(np.isin(y, (0.0, 1.0))):
        raise ValueError("detected must contain only binary outcomes 0 or 1.")
    if isinstance(min_replicates_per_level, (bool, np.bool_)) or not isinstance(
        min_replicates_per_level, (int, np.integer)
    ) or min_replicates_per_level < 2:
        raise ValueError("min_replicates_per_level must be an integer of at least 2.")
    if (
        isinstance(confidence_level, (bool, np.bool_))
        or not isinstance(confidence_level, Real)
        or not np.isfinite(confidence_level)
        or not 0.0 < confidence_level < 1.0
    ):
        raise ValueError("confidence_level must be a finite number strictly between 0 and 1.")

    levels, inverse, counts = np.unique(vaf, return_inverse=True, return_counts=True)
    if levels.size < 3:
        raise ValueError("at least 3 distinct VAF levels are required.")
    if np.any(counts < min_replicates_per_level):
        raise ValueError(
            f"each VAF level must have at least {min_replicates_per_level} replicate outcomes."
        )
    rates = np.asarray([y[inverse == index].mean() for index in range(levels.size)])
    if rates.min() >= 0.95 or rates.max() < 0.95:
        raise ValueError("empirical detection rates must bracket 0.95 across the VAF levels.")
    if np.unique(y).size < 2:
        raise ValueError("detected must include both 0 and 1 outcomes.")

    X = sm.add_constant(np.log10(vaf), has_constant="add")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", PerfectSeparationWarning)
            fit = sm.GLM(y, X, family=sm.families.Binomial(link=sm.families.links.Probit())).fit()
    except (PerfectSeparationError, PerfectSeparationWarning, np.linalg.LinAlgError) as exc:
        raise ValueError("probit fit is separated or unidentified; add replicated mixed outcomes near the transition.") from exc
    except (ValueError, FloatingPointError) as exc:
        raise ValueError(f"probit fit failed after input validation: {exc}") from exc

    params = np.asarray(fit.params, dtype=float)
    covariance = np.asarray(fit.cov_params(), dtype=float)
    if not fit.converged or params.shape != (2,) or covariance.shape != (2, 2) or not np.all(np.isfinite(covariance)):
        raise ValueError("probit fit did not converge to an identified two-parameter model.")
    intercept, slope = params
    if not np.isfinite(intercept) or not np.isfinite(slope) or slope <= 0.0:
        raise ValueError("probit fit must have a finite positive slope.")

    target_z = float(norm.ppf(0.95))
    log10_lod95 = float((target_z - intercept) / slope)
    gradient = np.asarray([-1.0 / slope, -(target_z - intercept) / slope**2])
    variance = float(gradient @ covariance @ gradient)
    if not np.isfinite(variance) or variance < 0.0:
        raise ValueError("probit LoD95 confidence interval is not identifiable from these data.")
    ci_z = float(norm.ppf(0.5 + confidence_level / 2.0))
    standard_error = float(np.sqrt(max(variance, 0.0)))
    ci_log10 = (log10_lod95 - ci_z * standard_error, log10_lod95 + ci_z * standard_error)
    lod95 = float(10**log10_lod95)
    ci = (float(10 ** ci_log10[0]), float(10 ** ci_log10[1]))
    if not np.isfinite(lod95) or not np.all(np.isfinite(ci)) or lod95 <= 0.0 or lod95 > 1.0:
        raise ValueError("fitted LoD95 is outside the physical VAF interval (0, 1].")

    return {
        "lod95_vaf": lod95,
        "confidence_level": float(confidence_level),
        "lod95_confidence_interval": [ci[0], ci[1]],
        "intercept": float(intercept),
        "slope": float(slope),
        "converged": bool(fit.converged),
        "n_observations": int(vaf.size),
        "n_levels": int(levels.size),
        "min_replicates_per_level": int(counts.min()),
        "deviance": float(fit.deviance),
        "pearson_chi2": float(fit.pearson_chi2),
        "residual_degrees_of_freedom": int(fit.df_resid),
        "empirical_levels": levels.tolist(),
        "empirical_detection_rates": rates.tolist(),
        "brackets_0_95": True,
    }


if __name__ == "__main__":
    rng = np.random.default_rng(42)
    blanks = np.clip(rng.normal(0.0002, 0.00005, size=20), 0.0, None)
    print(f"LoB (Gaussian shortcut, blank 95th pct) = {limit_of_blank(blanks):.6f} VAF")

    # Replicated, mixed outcomes bracket 95%; this is a deterministic teaching
    # fixture, not assay-validation evidence.
    levels = np.asarray([2e-4, 5e-4, 1e-3, 2e-3, 5e-3])
    positives = np.asarray([2, 8, 20, 34, 39])
    replicates = 40
    vaf_levels = np.repeat(levels, replicates)
    detected = np.concatenate(
        [np.r_[np.ones(count, dtype=int), np.zeros(replicates - count, dtype=int)] for count in positives]
    )
    result = lod95_probit(vaf_levels, detected)
    lower, upper = result["lod95_confidence_interval"]
    print("Replicated bracketing demo (contrived; not clinical validation evidence)")
    print(f"  LoD95 = {result['lod95_vaf']:.6f} VAF")
    print(f"  95% Wald CI = [{lower:.6f}, {upper:.6f}] VAF")
    print(
        f"  diagnostics: converged={result['converged']}, slope={result['slope']:.3f}, "
        f"n={result['n_observations']}, levels={result['n_levels']}, brackets_0.95={result['brackets_0_95']}"
    )
