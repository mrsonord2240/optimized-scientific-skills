#!/usr/bin/env python3
"""
Reference: PoPS (FinucaneLab/pops, git HEAD checked 2026-09-17); requires Python 3.9-3.11
with numpy/pandas/scipy/scikit-learn (`pip install -r pops/requirements.txt`); MAGMA
1.10+ gene-based output as input (see examples/magma_genebased.sh, Steps 1-2 -- Step 3
gene-set enrichment is NOT required for PoPS and should be skipped at locus scale).

Runs PoPS ridge (L2) regression on a MAGMA gene-based prefix + gene-feature matrix,
producing a per-gene polygenic priority score.

Windows-specific fix (verified 2026-09-17/18, see `references/pops.md`): MAGMA 1.10 on Windows writes `<prefix>.genes.out.txt` (extra `.txt`), but
`pops.py` hard-codes `<prefix>.genes.out` and raises FileNotFoundError if the exact
name is missing. This wrapper copies the file to the expected name first (no-op on
Linux/Mac, where MAGMA writes `.genes.out` directly).

Usage:
    python pops_run.py --pops_repo /path/to/pops --magma_prefix magma_gene \
        --gene_annot_path gene_annot.txt --feature_mat_prefix features \
        --control_features_path control.features --out_prefix pops_out
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def fix_windows_magma_output_name(magma_prefix: str) -> None:
    """Copy <prefix>.genes.out.txt to <prefix>.genes.out if the latter is missing."""
    expected = Path(f"{magma_prefix}.genes.out")
    windows_actual = Path(f"{magma_prefix}.genes.out.txt")
    if not expected.exists() and windows_actual.exists():
        shutil.copyfile(windows_actual, expected)
        print(f"[pops_run] Windows MAGMA output detected: copied "
              f"{windows_actual.name} -> {expected.name}")


def run_pops(pops_repo, magma_prefix, gene_annot_path, feature_mat_prefix,
             control_features_path, out_prefix, num_feature_chunks=1):
    fix_windows_magma_output_name(magma_prefix)
    pops_py = Path(pops_repo) / "pops.py"
    if not pops_py.exists():
        raise FileNotFoundError(f"pops.py not found under {pops_repo}")

    cmd = [
        sys.executable, str(pops_py),
        "--gene_annot_path", gene_annot_path,
        "--feature_mat_prefix", feature_mat_prefix,
        "--num_feature_chunks", str(num_feature_chunks),
        "--magma_prefix", magma_prefix,
        "--control_features_path", control_features_path,
        "--out_prefix", out_prefix,
    ]
    print("[pops_run] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    preds_path = f"{out_prefix}.preds"
    print(f"[pops_run] Done. Ranked genes from {preds_path}:")
    with open(preds_path) as f:
        header = f.readline().strip().split("\t")
        score_idx = header.index("PoPS_Score")
        rows = [line.strip().split("\t") for line in f]
    rows.sort(key=lambda r: float(r[score_idx]), reverse=True)
    for r in rows:
        print(f"  {r[0]:<20}{float(r[score_idx]):>14.6g}")

    print("\nCaveat: with too few genes / a single chromosome, PoPS's held-out-chromosome "
          "ridge CV has no fold to validate against and SELECTED_CV_ALPHA saturates, "
          "collapsing all scores toward 0 -- meaningful PoPS output needs genome-wide, "
          "multi-chromosome MAGMA input (see references/pops.md).")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pops_repo", required=True, help="Path to cloned FinucaneLab/pops repo")
    ap.add_argument("--magma_prefix", required=True)
    ap.add_argument("--gene_annot_path", required=True)
    ap.add_argument("--feature_mat_prefix", required=True)
    ap.add_argument("--control_features_path", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--num_feature_chunks", type=int, default=1)
    args = ap.parse_args()
    run_pops(args.pops_repo, args.magma_prefix, args.gene_annot_path,
              args.feature_mat_prefix, args.control_features_path,
              args.out_prefix, args.num_feature_chunks)
