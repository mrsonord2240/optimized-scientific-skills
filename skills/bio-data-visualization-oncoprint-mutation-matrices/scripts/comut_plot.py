#!/usr/bin/env python3
"""Build a cohort-complete comut.py plot from normalized long-format TSV files.

Inputs: mutation TSV (sample/category/value), cohort TSV (sample), and optional
clinical/TMB TSVs with sample/category/value. Usage: python scripts/comut_plot.py
--mutations mutations.tsv --cohort cohort.tsv --clinical clinical.tsv
--tmb tmb.tsv --output comut.pdf
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import pandas as pd
from comut import comut


ALTERATION_COLORS = {
    "Truncating": "#000000",
    "Missense": "#56B4E9",
    "Splice": "#CC79A7",
    "Amp": "#D55E00",
    "HomDel": "#0072B2",
    "Fusion": "#009E73",
}

DEFAULT_SAMPLE_LABEL_THRESHOLD = 50


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required column(s): {', '.join(missing)}")


def _normalize_sample_ids(
    values,
    label: str,
    *,
    require_unique: bool = False,
) -> list[str]:
    """Return stripped string IDs, rejecting missing or empty values."""
    normalized: list[str] = []
    for value in values:
        if pd.isna(value):
            raise ValueError(f"{label} sample IDs must not be empty or NA")
        sample = str(value).strip()
        if not sample:
            raise ValueError(f"{label} sample IDs must not be empty or NA")
        normalized.append(sample)
    if require_unique and len(normalized) != len(set(normalized)):
        raise ValueError(f"{label} sample IDs must be unique after normalization")
    return normalized


def _clinical_mapping(clinical: pd.DataFrame) -> dict[str, str]:
    values = sorted(clinical["value"].dropna().astype(str).unique())
    palette = plt.get_cmap("tab20")
    return {value: matplotlib.colors.to_hex(palette(i % palette.N)) for i, value in enumerate(values)}


def build_comut(
    mutations: pd.DataFrame,
    cohort_samples: list[str],
    clinical: pd.DataFrame | None = None,
    tmb: pd.DataFrame | None = None,
    top: int = 20,
    figsize: tuple[float, float] = (12, 8),
    sample_label_threshold: int = DEFAULT_SAMPLE_LABEL_THRESHOLD,
):
    """Return (CoMut object, descending-frequency genes, maximum TMB)."""
    if int(pd.__version__.split(".", maxsplit=1)[0]) >= 3:
        raise RuntimeError("comut 0.0.3 continuous tracks require pandas <3; use pandas 2.x.")
    _require_columns(mutations, {"sample", "category", "value"}, "mutations")
    if len(cohort_samples) == 0:
        raise ValueError("cohort sample IDs must not be empty")
    cohort_samples = _normalize_sample_ids(
        cohort_samples, "cohort", require_unique=True
    )
    if sample_label_threshold < 0:
        raise ValueError("sample_label_threshold must be zero or greater")
    cohort_set = set(cohort_samples)

    mutations = mutations.copy()
    mutations["sample"] = _normalize_sample_ids(mutations["sample"], "mutation")
    mutations[["category", "value"]] = mutations[["category", "value"]].astype(str)
    extra = sorted(set(mutations["sample"]).difference(cohort_set))
    if extra:
        raise ValueError(f"mutation samples absent from cohort: {', '.join(extra[:10])}")
    unknown = sorted(set(mutations["value"]).difference(ALTERATION_COLORS))
    if unknown:
        raise ValueError(f"unsupported alteration class(es): {', '.join(unknown)}")

    frequencies = mutations.groupby("category")["sample"].nunique().sort_values(ascending=False)
    top_genes = frequencies.index[:top].tolist()
    if not top_genes:
        raise ValueError("no mutation rows remain")
    mutations = mutations[mutations["category"].isin(top_genes)].drop_duplicates()

    burden = mutations.groupby("sample").size().reindex(cohort_samples, fill_value=0).astype(float)
    tmb_max = 0.0
    if tmb is not None:
        _require_columns(tmb, {"sample", "category", "value"}, "tmb")
        tmb = tmb.copy()
        tmb["sample"] = _normalize_sample_ids(tmb["sample"], "TMB")
        if set(tmb["sample"]).difference(cohort_set):
            raise ValueError("TMB data contains samples absent from the cohort")
        try:
            values = pd.to_numeric(tmb["value"], errors="raise")
        except (TypeError, ValueError) as exc:
            raise ValueError("TMB values must be finite and non-negative") from exc
        if values.isna().any() or any(
            not math.isfinite(float(value)) or float(value) < 0 for value in values
        ):
            raise ValueError("TMB values must be finite and non-negative")
        tmb_by_sample = pd.Series(values.to_numpy(), index=tmb["sample"]).groupby(level=0).max()
        burden = tmb_by_sample.reindex(cohort_samples, fill_value=0.0)
        tmb_max = float(burden.max())
        tmb = pd.DataFrame({
            "sample": cohort_samples,
            "category": "TMB",
            "value": burden.reindex(cohort_samples).to_numpy(dtype=float),
        })

    sample_order = sorted(cohort_samples, key=lambda sample: (-burden[sample], sample))
    plot = comut.CoMut()
    plot.samples = sample_order
    plot.add_categorical_data(
        data=mutations,
        name="Mutations",
        category_order=list(reversed(top_genes)),
        value_order=list(ALTERATION_COLORS),
        mapping=ALTERATION_COLORS,
    )

    if clinical is not None:
        _require_columns(clinical, {"sample", "category", "value"}, "clinical")
        clinical = clinical.copy()
        clinical["sample"] = _normalize_sample_ids(clinical["sample"], "clinical")
        if set(clinical["sample"]).difference(cohort_set):
            raise ValueError("clinical data contains samples absent from the cohort")
        plot.add_categorical_data(
            data=clinical,
            name="Clinical",
            mapping=_clinical_mapping(clinical),
        )

    if tmb is not None:
        plot.add_continuous_data(
            data=tmb,
            name="TMB",
            mapping="viridis",
            value_range=(0.0, max(1.0, tmb_max)),
        )

    plot.plot_comut(figsize=figsize)
    show_sample_labels = len(sample_order) <= sample_label_threshold
    for axis in plot.axes.values():
        axis.tick_params(axis="x", labelbottom=show_sample_labels)
    plot.figure.subplots_adjust(right=0.78)
    plot.add_unified_legend(
        axis_name="Mutations",
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
        frameon=False,
    )
    return plot, top_genes, tmb_max


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutations", type=Path, required=True)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--clinical", type=Path)
    parser.add_argument("--tmb", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--sample-label-threshold",
        type=int,
        default=DEFAULT_SAMPLE_LABEL_THRESHOLD,
        help=(
            "show sample IDs only when the cohort has at most this many samples "
            f"(default: {DEFAULT_SAMPLE_LABEL_THRESHOLD}; use 0 to hide all)"
        ),
    )
    args = parser.parse_args()

    cohort = _read_tsv(args.cohort)
    _require_columns(cohort, {"sample"}, "cohort")
    plot, _, _ = build_comut(
        mutations=_read_tsv(args.mutations),
        cohort_samples=cohort["sample"].tolist(),
        clinical=_read_tsv(args.clinical) if args.clinical else None,
        tmb=_read_tsv(args.tmb) if args.tmb else None,
        top=args.top,
        sample_label_threshold=args.sample_label_threshold,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot.figure.savefig(args.output, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
