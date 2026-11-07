"""Compute Holm/BH/Bonferroni-adjusted Mann-Whitney values and annotate them.

Inputs: CSV with group,value; output image; optional adjustment method.
Usage: python scripts/annotate_pairwise.py data.csv annotated.png holm
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator
from statsmodels.stats.multitest import multipletests


def annotate_pairwise(input_csv: str, output_plot: str, adjust_method: str = "holm") -> pd.DataFrame:
    df = pd.read_csv(input_csv)
    required = {"group", "value"}
    if not required.issubset(df.columns) or df[list(required)].isna().any().any():
        raise ValueError("input must contain non-missing group,value columns")
    groups = list(dict.fromkeys(df["group"].astype(str)))
    if len(groups) < 2:
        raise ValueError("at least two groups are required")
    pairs = list(itertools.combinations(groups, 2))
    test_results = [
        mannwhitneyu(
            df.loc[df["group"].astype(str) == left, "value"],
            df.loc[df["group"].astype(str) == right, "value"],
            alternative="two-sided",
            method="auto",
        )
        for left, right in pairs
    ]
    raw = [result.pvalue for result in test_results]
    effect_sizes = []
    for (left, right), result in zip(pairs, test_results, strict=True):
        n_left = int((df["group"].astype(str) == left).sum())
        n_right = int((df["group"].astype(str) == right).sum())
        effect_sizes.append(2 * float(result.statistic) / (n_left * n_right) - 1)
    method_map = {"BH": "fdr_bh", "bh": "fdr_bh", "fdr": "fdr_bh"}
    method = method_map.get(adjust_method, adjust_method)
    adjusted = multipletests(raw, method=method)[1]
    results = pd.DataFrame(
        {"group1": [p[0] for p in pairs], "group2": [p[1] for p in pairs],
         "test": "Mann-Whitney U", "p": raw, "p_adj": adjusted,
         "adjust_method": adjust_method, "family_size": len(pairs),
         "effect_type": "rank_biserial_r", "effect_size": effect_sizes}
    )

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sns.boxplot(data=df, x="group", y="value", order=groups, ax=ax)
    sns.stripplot(data=df, x="group", y="value", order=groups, color="black", alpha=0.45, ax=ax)
    annotator = Annotator(ax, pairs, data=df, x="group", y="value", order=groups)
    annotator.configure(test=None, text_format="star", line_height=0.02, text_offset=0.5)
    annotator.set_pvalues(adjusted)
    annotator.annotate()
    ax.set_title(f"Mann-Whitney; {adjust_method}-adjusted family ({len(pairs)} tests)")
    fig.tight_layout()
    fig.savefig(output_plot, dpi=160)
    plt.close(fig)

    results_path = Path(output_plot).with_suffix(".results.csv")
    results.to_csv(results_path, index=False)
    if not Path(output_plot).exists() or Path(output_plot).stat().st_size < 1000:
        raise RuntimeError("plot was not written")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv")
    parser.add_argument("output_plot")
    parser.add_argument("adjust_method", nargs="?", default="holm")
    args = parser.parse_args()
    annotate_pairwise(args.input_csv, args.output_plot, args.adjust_method)
