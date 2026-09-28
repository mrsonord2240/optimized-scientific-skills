"""Validate subject IDs, run a paired Wilcoxon test, and annotate a paired plot.

Inputs: CSV with subject_id,time,value; output image; optional adjustment method.
Usage: python scripts/annotate_paired.py paired.csv paired.png holm
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon
from statannotations.Annotator import Annotator
from statsmodels.stats.multitest import multipletests


def annotate_paired(input_csv: str, output_plot: str, adjust_method: str = "holm") -> pd.DataFrame:
    df = pd.read_csv(input_csv)
    required = ["subject_id", "time", "value"]
    if not set(required).issubset(df.columns) or df[required].isna().any().any():
        raise ValueError("input must contain non-missing subject_id,time,value columns")
    levels = list(dict.fromkeys(df["time"].astype(str)))
    if len(levels) != 2:
        raise ValueError("paired workflow requires exactly two time levels")
    if df.duplicated(["subject_id", "time"]).any():
        raise ValueError("duplicate subject/time rows")
    wide = df.pivot(index="subject_id", columns="time", values="value").reindex(columns=levels).sort_index()
    if wide.isna().any().any():
        raise ValueError("time levels do not contain identical complete subject IDs")

    raw_p = wilcoxon(wide[levels[0]], wide[levels[1]], alternative="two-sided", method="auto").pvalue
    method_map = {"BH": "fdr_bh", "bh": "fdr_bh", "fdr": "fdr_bh"}
    adjusted_p = multipletests([raw_p], method=method_map.get(adjust_method, adjust_method))[1][0]
    ordered = wide.reset_index().melt(id_vars="subject_id", value_vars=levels,
                                      var_name="time", value_name="value")

    fig, ax = plt.subplots(figsize=(5.5, 5))
    for _, subject in ordered.groupby("subject_id", sort=True):
        values = subject.set_index("time").reindex(levels)["value"]
        ax.plot(levels, values, color="#888888", alpha=0.45, linewidth=0.8)
    sns.boxplot(data=ordered, x="time", y="value", order=levels, ax=ax)
    annotator = Annotator(ax, [(levels[0], levels[1])], data=ordered,
                          x="time", y="value", order=levels)
    annotator.configure(test=None, text_format="simple")
    annotator.set_pvalues([adjusted_p])
    annotator.annotate()
    ax.set_title(f"Paired Wilcoxon; {adjust_method}-adjusted")
    fig.tight_layout()
    fig.savefig(output_plot, dpi=160)
    plt.close(fig)

    results = pd.DataFrame({"group1": [levels[0]], "group2": [levels[1]],
                            "p": [raw_p], "p_adj": [adjusted_p],
                            "adjust_method": [adjust_method], "n_pairs": [len(wide)]})
    results.to_csv(Path(output_plot).with_suffix(".results.csv"), index=False)
    if not Path(output_plot).exists() or Path(output_plot).stat().st_size < 1000:
        raise RuntimeError("plot was not written")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv")
    parser.add_argument("output_plot")
    parser.add_argument("adjust_method", nargs="?", default="holm")
    args = parser.parse_args()
    annotate_paired(args.input_csv, args.output_plot, args.adjust_method)
