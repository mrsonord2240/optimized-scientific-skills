#!/usr/bin/env python3
"""Validate, run, and summarize one RASQUAL invocation per cohort feature.

The runner targets the tested RASQUAL commit recorded in the skill reference.
It refuses partial/ambiguous manifests, checks native-double binary sizes, parses
the documented 25-column result contract, requires convergence, and adds a
Benjamini-Hochberg q-value across the actual rows emitted by the whole run.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable


OUTPUT_COLUMNS = [
    "feature_id", "rs_id", "chrom", "position", "ref", "alt", "allele_frequency",
    "hwe_chisq", "imputation_quality", "rasqual_log10_bh_q", "likelihood_ratio_chisq",
    "effect_pi", "error_delta", "mapping_bias_phi", "overdispersion", "snp_index",
    "n_feature_snps", "n_test_snps", "null_iterations", "alternative_iterations",
    "tie_location", "null_log_likelihood", "convergence_status", "fsnp_r2", "rsnp_r2",
]
MANIFEST_COLUMNS = [
    "feature_id", "region", "feature_index", "n_samples", "n_test_snps",
    "n_feature_snps", "exon_starts", "exon_ends", "lead_only",
]


class ContractError(RuntimeError):
    pass


def positive_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ContractError(f"{label} must be an integer") from exc
    if parsed < 1:
        raise ContractError(f"{label} must be positive")
    return parsed


def read_rectangular_table(path: Path, *, row_names: bool) -> tuple[int, int]:
    rows: list[list[str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row or all(not cell.strip() for cell in row):
                raise ContractError(f"{path}: blank row at line {line_number}")
            rows.append(row)
    if not rows:
        raise ContractError(f"{path}: no rows")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ContractError(f"{path}: ragged table")
    numeric_start = 1 if row_names else 0
    if width <= numeric_start:
        raise ContractError(f"{path}: no numeric values")
    for line_number, row in enumerate(rows, 1):
        if row_names and not row[0]:
            raise ContractError(f"{path}: empty feature ID at line {line_number}")
        try:
            values = [float(cell) for cell in row[numeric_start:]]
        except ValueError as exc:
            raise ContractError(f"{path}: non-numeric value at line {line_number}") from exc
        if not all(math.isfinite(value) for value in values):
            raise ContractError(f"{path}: non-finite value at line {line_number}")
    return len(rows), width - numeric_start


def prepare_binaries(args: argparse.Namespace) -> None:
    y = args.y.resolve(strict=True)
    k = args.k.resolve(strict=True)
    x = args.x.resolve(strict=True)
    y_features, y_samples = read_rectangular_table(y, row_names=True)
    k_features, k_samples = read_rectangular_table(k, row_names=True)
    x_samples, x_covariates = read_rectangular_table(x, row_names=False)
    if (y_features, y_samples) != (k_features, k_samples):
        raise ContractError("Y and K must have identical feature-by-sample dimensions")
    if x_samples != y_samples:
        raise ContractError("X rows must equal the Y/K sample count")
    if args.out.exists():
        raise ContractError(f"refusing to overwrite output directory: {args.out}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f"{args.out.name}.partial-", dir=args.out.parent))
    try:
        staged = {name: stage / f"{name}.txt" for name in ("Y", "K", "X")}
        for source, destination in zip((y, k, x), staged.values(), strict=True):
            shutil.copyfile(source, destination)
        converter = args.rasqual_source.resolve(strict=True) / "R" / "txt2bin.R"
        with converter.open("rb") as stdin:
            run = subprocess.run(
                [str(args.r), "--vanilla", "--quiet", "--args", *(str(path) for path in staged.values())],
                stdin=stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
        if run.returncode:
            raise ContractError(f"official txt2bin.R failed ({run.returncode}): {run.stderr.decode(errors='replace')}")
        expected = {"Y.bin": y_features * y_samples * 8, "K.bin": k_features * k_samples * 8,
                    "X.bin": x_samples * x_covariates * 8}
        for name, size in expected.items():
            actual = (stage / name).stat().st_size
            if actual != size:
                raise ContractError(f"{name}: expected {size} bytes, got {actual}")
        with (stage / "binary-contract.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["features", "samples", "covariates", "bytes_per_value"])
            writer.writerow([y_features, y_samples, x_covariates, 8])
        stage.rename(args.out)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise


def read_manifest(path: Path, expected_features: int) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != MANIFEST_COLUMNS:
            raise ContractError(f"manifest columns must be exactly: {','.join(MANIFEST_COLUMNS)}")
        rows = list(reader)
    if len(rows) != expected_features:
        raise ContractError(f"manifest must contain exactly {expected_features} feature rows")
    indices = [positive_int(row["feature_index"], "feature_index") for row in rows]
    if sorted(indices) != list(range(1, expected_features + 1)) or len(set(indices)) != len(indices):
        raise ContractError("feature_index must cover 1..feature-count exactly once")
    ids = [row["feature_id"] for row in rows]
    if any(not feature_id for feature_id in ids) or len(set(ids)) != len(ids):
        raise ContractError("feature_id must be non-empty and unique")
    for row in rows:
        for key in ("n_samples", "n_test_snps", "n_feature_snps"):
            positive_int(row[key], key)
        if row["lead_only"] not in {"true", "false"}:
            raise ContractError("lead_only must be true or false")
        if ":" not in row["region"] or "-" not in row["region"]:
            raise ContractError(f"invalid tabix region: {row['region']}")
        for key in ("exon_starts", "exon_ends"):
            values = [positive_int(value, key) for value in row[key].split(",")]
            if values != sorted(values):
                raise ContractError(f"{key} must be ascending")
        if len(row["exon_starts"].split(",")) != len(row["exon_ends"].split(",")):
            raise ContractError("exon_starts and exon_ends lengths differ")
    return sorted(rows, key=lambda row: int(row["feature_index"]))


def assert_binary_size(path: Path, expected: int, label: str) -> None:
    actual = path.resolve(strict=True).stat().st_size
    if actual != expected:
        raise ContractError(f"{label}: expected {expected} bytes, got {actual}")


def chisq1_survival(chisq: float) -> float:
    if not math.isfinite(chisq) or chisq < 0:
        raise ContractError(f"invalid likelihood-ratio chi-square: {chisq}")
    return math.erfc(math.sqrt(chisq / 2.0))


def benjamini_hochberg(p_values: Iterable[float]) -> list[float]:
    values = list(p_values)
    if not values:
        return []
    if any(not math.isfinite(value) or value < 0 or value > 1 for value in values):
        raise ContractError("p-values must be finite within [0,1]")
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [1.0] * len(values)
    running = 1.0
    count = len(values)
    for rank_index in range(count - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[index] * count / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def parse_output(text: str, expected_feature: str) -> list[list[str]]:
    rows = list(csv.reader((line for line in text.splitlines() if line.strip()), delimiter="\t"))
    if not rows:
        raise ContractError(f"{expected_feature}: RASQUAL emitted no rows")
    for row in rows:
        if len(row) != 25:
            raise ContractError(f"{expected_feature}: expected 25 columns, got {len(row)}")
        if row[0] != expected_feature:
            raise ContractError(f"{expected_feature}: output feature is {row[0]}")
        try:
            phi = float(row[13])
            convergence = int(row[22])
        except ValueError as exc:
            raise ContractError(f"{expected_feature}: invalid phi or convergence field") from exc
        if not math.isfinite(phi) or not 0 <= phi <= 1:
            raise ContractError(f"{expected_feature}: phi outside [0,1]")
        if convergence != 0:
            raise ContractError(f"{expected_feature}: convergence status is {convergence}")
    return rows


def run_cohort(args: argparse.Namespace) -> None:
    manifest = read_manifest(args.manifest.resolve(strict=True), args.feature_count)
    n_samples = {positive_int(row["n_samples"], "n_samples") for row in manifest}
    if len(n_samples) != 1:
        raise ContractError("all manifest rows must use the same n_samples")
    samples = n_samples.pop()
    assert_binary_size(args.y, args.feature_count * samples * 8, "Y.bin")
    assert_binary_size(args.k, args.feature_count * samples * 8, "K.bin")
    assert_binary_size(args.x, samples * args.covariates * 8, "X.bin")
    for executable in (args.tabix, args.rasqual):
        if not Path(executable).resolve(strict=True).is_file():
            raise ContractError(f"missing executable: {executable}")
    args.vcf.resolve(strict=True)
    if args.output.exists():
        raise ContractError(f"refusing to overwrite output: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    all_rows: list[list[str]] = []
    for row in manifest:
        tabix = subprocess.run(
            [str(args.tabix), str(args.vcf), row["region"]],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        if tabix.returncode:
            raise ContractError(f"{row['feature_id']}: tabix failed: {tabix.stderr.decode(errors='replace')}")
        command = [
            str(args.rasqual), "-y", str(args.y), "-k", str(args.k), "-x", str(args.x),
            "-n", row["n_samples"], "-j", row["feature_index"], "-l", row["n_test_snps"],
            "-m", row["n_feature_snps"], "-s", row["exon_starts"], "-e", row["exon_ends"],
            "-f", row["feature_id"], "-z",
        ]
        if row["lead_only"] == "true":
            command.append("-t")
        result = subprocess.run(command, input=tabix.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if result.returncode:
            raise ContractError(f"{row['feature_id']}: RASQUAL failed: {result.stderr.decode(errors='replace')}")
        all_rows.extend(parse_output(result.stdout.decode(), row["feature_id"]))

    raw_p = [chisq1_survival(float(row[10])) for row in all_rows]
    family_q = benjamini_hochberg(raw_p)
    fd, temp_name = tempfile.mkstemp(prefix=f".{args.output.name}.", dir=args.output.parent, text=True)
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([*OUTPUT_COLUMNS, "cohort_raw_p", "cohort_bh_q"])
            for row, p_value, q_value in zip(all_rows, raw_p, family_q, strict=True):
                writer.writerow([*row, f"{p_value:.17g}", f"{q_value:.17g}"])
        Path(temp_name).replace(args.output)
    except Exception:
        Path(temp_name).unlink(missing_ok=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)
    prepare = sub.add_parser("prepare", help="validate text matrices and run official txt2bin.R")
    prepare.add_argument("--rasqual-source", type=Path, required=True)
    prepare.add_argument("--r", type=Path, required=True, help="R executable, not Rscript")
    prepare.add_argument("--y", type=Path, required=True)
    prepare.add_argument("--k", type=Path, required=True)
    prepare.add_argument("--x", type=Path, required=True)
    prepare.add_argument("--out", type=Path, required=True)
    prepare.set_defaults(func=prepare_binaries)

    run = sub.add_parser("run", help="execute each manifest feature once and summarize the cohort family")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--feature-count", type=int, required=True)
    run.add_argument("--covariates", type=int, required=True)
    run.add_argument("--y", type=Path, required=True)
    run.add_argument("--k", type=Path, required=True)
    run.add_argument("--x", type=Path, required=True)
    run.add_argument("--vcf", type=Path, required=True)
    run.add_argument("--tabix", type=Path, required=True)
    run.add_argument("--rasqual", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    run.set_defaults(func=run_cohort)
    return parser


def main() -> int:
    try:
        args = build_parser().parse_args()
        if getattr(args, "feature_count", 1) < 1 or getattr(args, "covariates", 1) < 1:
            raise ContractError("feature-count and covariates must be positive")
        args.func(args)
        return 0
    except ContractError as exc:
        print(f"contract error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
