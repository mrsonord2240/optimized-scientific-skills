"""Validated limit-of-blank and replicated LoD95 probit calculations."""

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


def main():
    import argparse
    import csv
    import sys

    import pandas as pd

    ap = argparse.ArgumentParser(description="LoD95 probit fit on a replicated dilution series, optional LoB from blanks.")
    ap.add_argument("dilution_csv", help="replicate rows (vaf, detected 0/1) or summary rows (vaf, n_replicates, n_detected)")
    ap.add_argument("--blanks", help="CSV with a 'signal' column of blank-replicate VAF-scale readings")
    ap.add_argument("--vaf-col", default="vaf")
    ap.add_argument("--detected-col", default="detected")
    ap.add_argument("--confidence", type=float, default=0.95)
    ap.add_argument("--min-replicates", type=int, default=3, help="replicates required at every VAF level")
    ap.add_argument("--out", default="lod95_results.tsv")
    a = ap.parse_args()
    df = pd.read_csv(a.dilution_csv)
    if "n_replicates" in df.columns and "n_detected" in df.columns:
        if (df["n_detected"] > df["n_replicates"]).any() or (df["n_detected"] < 0).any():
            sys.exit("REFUSED: n_detected must be between 0 and n_replicates")
        vaf, det = [], []
        for _, r in df.iterrows():
            k, n = int(r["n_detected"]), int(r["n_replicates"])
            vaf += [r[a.vaf_col]] * n
            det += [1] * k + [0] * (n - k)
    elif a.vaf_col in df.columns and a.detected_col in df.columns:
        vaf, det = df[a.vaf_col].tolist(), df[a.detected_col].tolist()
    else:
        sys.exit(f"REFUSED: need columns {a.vaf_col!r} and {a.detected_col!r}, or n_replicates and n_detected; found {list(df.columns)}")
    if max(vaf) > 1:
        sys.exit("REFUSED: vaf is a fraction, not a percent (0.1% is 0.001); the largest value here exceeds 1.")
    try:
        res = lod95_probit(vaf, det, confidence_level=a.confidence, min_replicates_per_level=a.min_replicates)
        rows = [("lod95_vaf", res["lod95_vaf"]),
                ("lod95_ci_low", res["lod95_confidence_interval"][0]), ("lod95_ci_high", res["lod95_confidence_interval"][1]),
                ("confidence_level", res["confidence_level"]), ("slope", res["slope"]), ("intercept", res["intercept"]),
                ("converged", res["converged"]), ("n_observations", res["n_observations"]), ("n_levels", res["n_levels"]),
                ("min_replicates_per_level", res["min_replicates_per_level"]), ("deviance", res["deviance"]),
                ("pearson_chi2", res["pearson_chi2"]), ("residual_df", res["residual_degrees_of_freedom"])]
        rows += [(f"detection_rate_at_{lv:g}", rt) for lv, rt in zip(res["empirical_levels"], res["empirical_detection_rates"])]
        if a.blanks:
            b = pd.read_csv(a.blanks)
            if "signal" not in b.columns:
                sys.exit("REFUSED: blanks CSV needs a 'signal' column")
            rows.append(("lob_gaussian", limit_of_blank(b["signal"].tolist())))
    except ValueError as exc:
        sys.exit(f"REFUSED: {exc}")
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerows(rows)
    for k, v in rows:
        print(f"{k}\t{v}")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
