"""Regression test for scripts/comut_plot.py using the audit's planted cohort."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import subprocess
import sys

import pandas as pd


CLASS_MAP = {
    "Missense_Mutation": "Missense",
    "In_Frame_Ins": "Missense",
    "In_Frame_Del": "Missense",
    "Nonsense_Mutation": "Truncating",
    "Frame_Shift_Ins": "Truncating",
    "Frame_Shift_Del": "Truncating",
    "Nonstop_Mutation": "Truncating",
    "Translation_Start_Site": "Truncating",
    "Splice_Site": "Splice",
}


def expect_value_error(label: str, call) -> None:
    try:
        call()
    except ValueError:
        return
    raise AssertionError(f"{label} was accepted")


def rendered_texts(text_artists, figure) -> list[str]:
    """Return visible, non-empty text whose rendered bounding box has area."""
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    return [
        artist.get_text()
        for artist in text_artists
        if artist.get_visible()
        and artist.get_text()
        and artist.get_window_extent(renderer).width > 0
        and artist.get_window_extent(renderer).height > 0
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("maf", type=Path)
    parser.add_argument("clinical", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    script = Path(__file__).parents[1] / "scripts" / "comut_plot.py"
    spec = importlib.util.spec_from_file_location("comut_plot", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    maf = pd.read_csv(args.maf, sep="\t")
    clinical = pd.read_csv(args.clinical, sep="\t")
    mapped = maf["Variant_Classification"].map(CLASS_MAP)
    mutations = pd.DataFrame({
        "sample": maf.loc[mapped.notna(), "Tumor_Sample_Barcode"],
        "category": maf.loc[mapped.notna(), "Hugo_Symbol"],
        "value": mapped.dropna(),
    }).drop_duplicates()
    cohort = clinical["Tumor_Sample_Barcode"].astype(str).tolist()
    clinical_long = clinical.rename(columns={
        "Tumor_Sample_Barcode": "sample", "Subtype": "value"
    }).assign(category="Subtype")[["sample", "category", "value"]]
    tmb_values = maf.groupby("Tumor_Sample_Barcode").size().reindex(cohort, fill_value=0)
    tmb = pd.DataFrame({"sample": cohort, "category": "TMB", "value": tmb_values.to_numpy()})

    args.output_dir.mkdir(parents=True, exist_ok=True)
    mutation_path = args.output_dir / "mutations.tsv"
    cohort_path = args.output_dir / "cohort.tsv"
    clinical_path = args.output_dir / "clinical.tsv"
    tmb_path = args.output_dir / "tmb.tsv"
    cli_output = args.output_dir / "comut-cli.pdf"
    mutations.to_csv(mutation_path, sep="\t", index=False)
    pd.DataFrame({"sample": cohort}).to_csv(cohort_path, sep="\t", index=False)
    clinical_long.to_csv(clinical_path, sep="\t", index=False)
    tmb.to_csv(tmb_path, sep="\t", index=False)
    subprocess.run([
        sys.executable, str(script),
        "--mutations", str(mutation_path),
        "--cohort", str(cohort_path),
        "--clinical", str(clinical_path),
        "--tmb", str(tmb_path),
        "--output", str(cli_output),
        "--top", "6",
    ], check=True)
    assert cli_output.stat().st_size > 1_000

    plot, top_genes, tmb_max = module.build_comut(
        mutations, cohort, clinical=clinical_long, tmb=tmb, top=6
    )
    expected_samples = sorted(cohort, key=lambda sample: (-tmb_values[sample], sample))
    ticks = [tick.get_text() for tick in plot.axes["Mutations"].get_yticklabels()]
    assert plot.samples == expected_samples
    assert len(plot.samples) == len(cohort)
    assert ticks[-1] == top_genes[0]
    assert tmb_max == float(tmb_values.max())
    legend = plot.axes["Mutations"].get_legend()
    assert legend is not None
    legend_text = rendered_texts(legend.get_texts(), plot.figure)
    observed_classes = set(mutations["value"])
    observed_clinical = set(clinical_long["value"].dropna().astype(str))
    assert {"Mutations", "Clinical"}.issubset(legend_text)
    assert observed_classes.issubset(legend_text)
    assert observed_clinical.issubset(legend_text)
    small_sample_labels = rendered_texts(
        plot.axes["Mutations"].get_xticklabels(), plot.figure
    )
    assert small_sample_labels == expected_samples
    direct_output = args.output_dir / "comut-direct.png"
    plot.figure.savefig(direct_output, dpi=100, bbox_inches="tight")
    assert direct_output.stat().st_size > 1_000

    dense_cohort = cohort + [f"EMPTY_{index:03d}" for index in range(31)]
    dense_plot, _, _ = module.build_comut(mutations, dense_cohort, top=6)
    assert not rendered_texts(
        dense_plot.axes["Mutations"].get_xticklabels(), dense_plot.figure
    )
    dense_legend = dense_plot.axes["Mutations"].get_legend()
    assert dense_legend is not None
    assert observed_classes.issubset(rendered_texts(dense_legend.get_texts(), dense_plot.figure))
    dense_output = args.output_dir / "comut-dense.png"
    dense_plot.figure.savefig(dense_output, dpi=100, bbox_inches="tight")
    assert dense_output.stat().st_size > 1_000

    for bad_id in ["", "   ", None, pd.NA, float("nan")]:
        expect_value_error(
            f"cohort ID {bad_id!r}",
            lambda bad_id=bad_id: module.build_comut(mutations, cohort + [bad_id]),
        )

    for bad_tmb in ["", pd.NA, float("inf"), float("-inf"), -1]:
        bad = tmb.copy()
        bad["value"] = bad["value"].astype(object)
        bad.loc[0, "value"] = bad_tmb
        expect_value_error(
            f"TMB value {bad_tmb!r}",
            lambda bad=bad: module.build_comut(mutations, cohort, tmb=bad, top=6),
        )

    print(
        "PASS: cohort/order/TMB range, rendered legends, dense label policy, "
        "and invalid ID/TMB rejection"
    )


if __name__ == "__main__":
    main()
