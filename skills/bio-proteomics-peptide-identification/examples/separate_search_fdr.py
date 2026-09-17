"""FDR for SEPARATE target and decoy searches: pi0 * D / T (Kall et al. 2008).

The concatenated-competition estimator (decoys + 1) / targets in SKILL.md is
WRONG here. In a separate search every spectrum gets both a target hit and a
decoy hit, so no competition resolved which one wins, and the decoy count
estimates the number of INCORRECT target PSMs directly -- scaled by pi0, the
proportion of target PSMs that are incorrect.

Inputs are two tables with one best hit per spectrum each, both with columns
'scan' and 'score' (higher is better; use -log10(E-value) for E-value scores).

Run:  python separate_search_fdr.py target.tsv decoy.tsv
Checked on Python 3.12, pandas 2.2+, numpy 1.26+.
"""

import sys

import numpy as np
import pandas as pd


def estimate_pi0(target_scores, decoy_scores):
    """Kall et al. 2008 (JPR 7:29-34) median-decoy estimator of pi0.

    Half of the incorrect target PSMs are expected to score below the median
    decoy score, so pi0 = 2 * (fraction of targets below it). Capped at 1,
    which is the conservative fallback when the estimate runs high.
    """
    median_decoy = np.median(decoy_scores)
    below = float(np.mean(np.asarray(target_scores) < median_decoy))
    return min(1.0, 2.0 * below)


def separate_search_qvalues(targets, decoys, pi0=None):
    """q-values for a separate target/decoy search.

    targets, decoys: DataFrames with 'scan' and 'score', one row per spectrum.
    pi0=None estimates pi0 from the score distributions; pass pi0=1.0 for the
    valid-but-conservative bound.
    """
    for name, df in (("targets", targets), ("decoys", decoys)):
        if df["scan"].duplicated().any():
            raise ValueError(f"{name}: more than one row per scan -- keep the best hit per spectrum first")
    if len(decoys) == 0:
        raise ValueError("no decoy PSMs: a separate-search FDR needs the decoy table, not a decoy-filtered one")

    if pi0 is None:
        pi0 = estimate_pi0(targets["score"].to_numpy(), decoys["score"].to_numpy())

    t = targets.sort_values("score", ascending=False).reset_index(drop=True)
    decoy_sorted = np.sort(decoys["score"].to_numpy())
    # decoys at or above each target score, without an O(n^2) scan
    n_decoy_above = len(decoy_sorted) - np.searchsorted(decoy_sorted, t["score"].to_numpy(), side="left")
    n_target_above = np.arange(1, len(t) + 1)

    t["fdr"] = np.minimum(1.0, pi0 * n_decoy_above / n_target_above)
    t["qvalue"] = t["fdr"][::-1].cummin()[::-1]
    return t, pi0


def main(target_path, decoy_path):
    targets = pd.read_csv(target_path, sep="\t")
    decoys = pd.read_csv(decoy_path, sep="\t")
    # these synthetic tables carry one row per scan already; a real export needs
    # targets = targets.sort_values('score', ascending=False).drop_duplicates('scan')
    scored, pi0 = separate_search_qvalues(targets, decoys)
    kept = scored[scored["qvalue"] <= 0.01]
    print(f"targets {len(targets)}  decoys {len(decoys)}  pi0-hat {pi0:.3f}")
    print(f"kept at q <= 0.01: {len(kept)}")

    conservative, _ = separate_search_qvalues(targets, decoys, pi0=1.0)
    print(f"kept at q <= 0.01 with pi0 = 1 (conservative): {int((conservative['qvalue'] <= 0.01).sum())}")

    if "is_correct" in scored.columns:  # synthetic data with ground truth
        print(f"true FDP of the pi0-hat list: {1 - kept['is_correct'].mean():.4f}")
        cons_kept = conservative[conservative["qvalue"] <= 0.01]
        print(f"true FDP of the pi0 = 1 list: {1 - cons_kept['is_correct'].mean():.4f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python separate_search_fdr.py target.tsv decoy.tsv")
    main(sys.argv[1], sys.argv[2])
