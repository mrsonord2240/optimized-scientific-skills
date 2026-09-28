"""Focused checks for adjusted Python labels and subject-ID-safe pairing."""

from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests


SKILL_DIR = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


pairwise_mod = load_module("annotate_pairwise", SKILL_DIR / "scripts" / "annotate_pairwise.py")
paired_mod = load_module("annotate_paired", SKILL_DIR / "scripts" / "annotate_paired.py")

with tempfile.TemporaryDirectory(prefix="statanno-regression-") as tmp:
    scratch = Path(tmp)
    data_dir = SKILL_DIR / "examples" / "data"
    results = pairwise_mod.annotate_pairwise(
        str(data_dir / "three_group.csv"), str(scratch / "pairwise.png"), "holm"
    )
    expected = multipletests(results["p"], method="holm")[1]
    assert max(abs(results["p_adj"] - expected)) < 1e-12

    paired = paired_mod.annotate_paired(
        str(data_dir / "paired.csv"), str(scratch / "paired.png"), "holm"
    )
    assert paired.loc[0, "p_adj"] < 0.01

    incomplete = pd.read_csv(data_dir / "paired.csv").iloc[1:]
    incomplete_path = scratch / "paired-incomplete.csv"
    incomplete.to_csv(incomplete_path, index=False)
    try:
        paired_mod.annotate_paired(str(incomplete_path), str(scratch / "should-not-exist.png"))
    except ValueError:
        pass
    else:
        raise AssertionError("incomplete subject pairs were not rejected")

print("Python regression checks passed")
