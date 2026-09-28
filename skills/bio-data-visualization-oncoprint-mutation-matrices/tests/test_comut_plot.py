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
    direct_output = args.output_dir / "comut-direct.png"
    plot.figure.savefig(direct_output, dpi=100, bbox_inches="tight")
    assert direct_output.stat().st_size > 1_000
    print("PASS: full cohort, burden order, top gene placement, and data-derived TMB range")


if __name__ == "__main__":
    main()
