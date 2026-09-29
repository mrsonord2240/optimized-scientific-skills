#!/usr/bin/env python3
"""Run the pinned RASQUAL CLI once per validated feature row.

The driver owns the 1-based ``-j`` matrix row and writes a feature-level BH
summary from the maximum likelihood-ratio chi-square reported for each feature.
It invokes tabix and RASQUAL without a shell so paths and feature names are not
reinterpreted.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Feature:
    index: int
    name: str
    chrom: str
    start: int
    end: int
    testing_snps: int
    feature_snps: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--offsets", required=True, type=Path)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=int)
    parser.add_argument("--matrix-rows", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--rasqual", default="rasqual")
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--cis-window", default=500_000, type=int)
    parser.add_argument(
        "--chi-square-column", default=11, type=int,
        help="One-based RASQUAL output column containing the joint chi-square (default: 11)",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def fail(message: str) -> "NoReturn":
    raise ValueError(message)


def read_features(path: Path, matrix_rows: int) -> tuple[list[Feature], list[str]]:
    features: list[Feature] = []
    skipped: list[str] = []
    names: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) != 6:
                fail(f"features line {line_number} must have exactly six tab-delimited fields")
            name, chrom = row[0], row[1]
            if not name or name in names:
                fail(f"features line {line_number} has an empty or duplicate name: {name!r}")
            names.add(name)
            try:
                start, end, n_testing, n_feature = map(int, row[2:])
            except ValueError as exc:
                fail(f"features line {line_number} has non-integer coordinates or SNP counts: {exc}")
            if start < 1 or end < start:
                fail(f"features line {line_number} requires 1 <= start <= end")
            if n_testing < 0 or n_feature < 0 or n_feature > n_testing:
                fail(f"features line {line_number} has inconsistent SNP counts")
            index = len(features) + len(skipped) + 1
            if index > matrix_rows:
                fail(f"feature row {index} exceeds declared matrix rows {matrix_rows}")
            if n_feature == 0:
                skipped.append(name)
                continue
            features.append(Feature(index, name, chrom, start, end, n_testing, n_feature))
    if not features:
        fail("no executable feature has at least one feature SNP")
    return features, skipped


def command_for(feature: Feature, args: argparse.Namespace) -> tuple[list[str], list[str]]:
    region_start = max(1, feature.start - args.cis_window)
    region_end = feature.end + args.cis_window
    tabix = [args.tabix, str(args.vcf), f"{feature.chrom}:{region_start}-{region_end}"]
    rasqual = [
        args.rasqual, "-y", str(args.counts), "-k", str(args.offsets),
        "-n", str(args.samples), "-j", str(feature.index),
        "-l", str(feature.testing_snps), "-m", str(feature.feature_snps),
        "-s", str(feature.start), "-e", str(feature.end), "-f", feature.name,
    ]
    return tabix, rasqual


def validate_runtime(args: argparse.Namespace) -> None:
    if args.samples < 1 or args.matrix_rows < 1 or args.cis_window < 0:
        fail("samples and matrix-rows must be positive; cis-window cannot be negative")
    if args.chi_square_column < 1:
        fail("chi-square-column is one-based and must be positive")
    for path in (args.features, args.counts, args.offsets, args.vcf):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"required non-empty input is missing: {path}")
    expected_matrix_bytes = args.samples * args.matrix_rows * 8
    for label, path in (("counts", args.counts), ("offsets", args.offsets)):
        if path.stat().st_size != expected_matrix_bytes:
            fail(
                f"{label} matrix has {path.stat().st_size} bytes; expected exactly "
                f"{args.samples} samples x {args.matrix_rows} rows x 8 bytes = "
                f"{expected_matrix_bytes}"
            )
    opener = gzip.open if args.vcf.suffix == ".gz" else open
    try:
        with opener(args.vcf, mode="rt", encoding="utf-8") as handle:
            header = next((line for line in handle if line.startswith("#CHROM\t")), None)
    except (EOFError, gzip.BadGzipFile, UnicodeDecodeError) as exc:
        fail(f"could not read VCF header: {exc}")
    if header is None:
        fail("VCF has no #CHROM header")
    vcf_samples = max(0, len(header.rstrip("\r\n").split("\t")) - 9)
    if vcf_samples != args.samples:
        fail(f"VCF has {vcf_samples} samples but --samples declares {args.samples}")
    if not args.dry_run:
        for executable in (args.rasqual, args.tabix):
            if shutil.which(executable) is None:
                fail(f"required executable is unavailable: {executable}")
        if not Path(f"{args.vcf}.tbi").is_file() and not Path(f"{args.vcf}.csi").is_file():
            fail(f"indexed VCF is required (.tbi or .csi): {args.vcf}")
        if args.output_dir.exists():
            fail(f"output directory already exists: {args.output_dir}")


def run_feature(feature: Feature, args: argparse.Namespace, output: Path) -> list[float]:
    tabix_cmd, rasqual_cmd = command_for(feature, args)
    with output.open("wb") as result:
        tabix = subprocess.Popen(tabix_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert tabix.stdout is not None
        rasqual = subprocess.Popen(
            rasqual_cmd, stdin=tabix.stdout, stdout=result, stderr=subprocess.PIPE
        )
        tabix.stdout.close()
        rasqual_stderr = rasqual.communicate()[1]
        tabix_stderr = tabix.communicate()[1]
    if tabix.returncode != 0:
        fail(f"tabix failed for {feature.name}: {tabix_stderr.decode(errors='replace').strip()}")
    if rasqual.returncode != 0:
        fail(f"RASQUAL failed for {feature.name}: {rasqual_stderr.decode(errors='replace').strip()}")
    chi_squares: list[float] = []
    with output.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) < args.chi_square_column:
                fail(f"RASQUAL row {line_number} for {feature.name} lacks chi-square column")
            try:
                value = float(fields[args.chi_square_column - 1])
            except ValueError as exc:
                fail(f"invalid chi-square for {feature.name} row {line_number}: {exc}")
            if not math.isfinite(value) or value < 0:
                fail(f"non-finite or negative chi-square for {feature.name} row {line_number}")
            chi_squares.append(value)
    if not chi_squares:
        fail(f"RASQUAL returned no rows for {feature.name}")
    return chi_squares


def bh_adjust(p_values: list[float]) -> list[float]:
    count = len(p_values)
    adjusted = [1.0] * count
    running = 1.0
    for rank, index in reversed(list(enumerate(sorted(range(count), key=p_values.__getitem__), start=1))):
        running = min(running, p_values[index] * count / rank)
        adjusted[index] = running
    return adjusted


def main() -> int:
    args = parse_args()
    try:
        validate_runtime(args)
        features, skipped = read_features(args.features, args.matrix_rows)
        if args.dry_run:
            print(json.dumps({
                "features": [
                    {
                        "index": feature.index,
                        "name": feature.name,
                        "tabix": command_for(feature, args)[0],
                        "rasqual": command_for(feature, args)[1],
                    }
                    for feature in features
                ],
                "skipped_zero_feature_snps": skipped,
                "multiplicity": "BH across every reported feature-variant joint test",
            }, indent=2))
            return 0
        args.output_dir.mkdir(parents=True)
        rows = []
        for feature in features:
            chi_squares = run_feature(
                feature, args, args.output_dir / f"{feature.index:06d}.rasqual.tsv"
            )
            for rasqual_row, chi_square in enumerate(chi_squares, start=1):
                p_value = math.erfc(math.sqrt(chi_square / 2.0))
                rows.append([feature.index, feature.name, rasqual_row, chi_square, p_value])
        adjusted = bh_adjust([row[4] for row in rows])
        with (args.output_dir / "association_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "feature_index", "feature_name", "rasqual_row", "joint_chi_square",
                "p_value", "adj_p",
            ])
            for row, adj_p in zip(rows, adjusted):
                writer.writerow([*row, adj_p])
        with (args.output_dir / "skipped_zero_feature_snps.txt").open("w", encoding="utf-8") as handle:
            handle.write("\n".join(skipped) + ("\n" if skipped else ""))
        return 0
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"RASQUAL driver error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
