#!/usr/bin/env python3
"""Validate current Hyb rows and retain only assignments stable across runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

MIRNA_RE = re.compile(r"(?:^|_)microRNA(?:_|$)", re.IGNORECASE)
TARGET_RE = re.compile(r"(?:^|_)mRNA(?:_|$)", re.IGNORECASE)
SITE_HEADER = [
    "read_id", "sequence", "folding_energy", "mirna_id", "target_id",
    "mirna_read_start", "mirna_read_end", "mirna_reference_start",
    "mirna_reference_end", "mirna_alignment_score", "target_read_start",
    "target_read_end", "target_reference_start", "target_reference_end",
    "target_alignment_score", "source_orientation", "supporting_runs",
    "expression_value", "expression_unit", "expression_source",
]


@dataclass(frozen=True)
class Segment:
    name: str
    read_start: str
    read_end: str
    reference_start: str
    reference_end: str
    score: str


@dataclass(frozen=True)
class Site:
    read_id: str
    sequence: str
    energy: str
    mirna: Segment
    target: Segment
    orientation: str

    def assignment_key(self) -> tuple[str, ...]:
        return (
            self.sequence, self.energy, self.mirna.name, self.target.name,
            self.mirna.read_start, self.mirna.read_end,
            self.mirna.reference_start, self.mirna.reference_end,
            self.mirna.score, self.target.read_start, self.target.read_end,
            self.target.reference_start, self.target.reference_end, self.target.score,
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def segment(row: list[str], start: int) -> Segment:
    return Segment(*row[start : start + 6])


def parse_site(row: list[str], path: Path, line_no: int) -> Site:
    if len(row) != 16:
        raise ValueError(f"{path}:{line_no}: expected 16 tab fields, found {len(row)}")
    if not row[0] or not row[1] or not row[3] or not row[9]:
        raise ValueError(f"{path}:{line_no}: required Hyb fields are empty")
    first, second = segment(row, 3), segment(row, 9)
    first_mirna = bool(MIRNA_RE.search(first.name)); second_mirna = bool(MIRNA_RE.search(second.name))
    first_target = bool(TARGET_RE.search(first.name)); second_target = bool(TARGET_RE.search(second.name))
    if first_mirna and second_target:
        mirna, target, orientation = first, second, "mirna-first"
    elif second_mirna and first_target:
        mirna, target, orientation = second, first, "target-first"
    else:
        raise ValueError(
            f"{path}:{line_no}: fields 4/10 are not one microRNA and one mRNA: "
            f"{first.name!r}, {second.name!r}"
        )
    for label, current in (("mirna", mirna), ("target", target)):
        values = (current.read_start, current.read_end, current.reference_start, current.reference_end)
        if any(not value.isdigit() or int(value) < 1 for value in values):
            raise ValueError(f"{path}:{line_no}: invalid {label} coordinate fields: {values!r}")
        if int(current.read_start) > int(current.read_end) or int(current.reference_start) > int(current.reference_end):
            raise ValueError(f"{path}:{line_no}: reversed {label} coordinate interval")
        try:
            float(current.score)
        except ValueError as exc:
            raise ValueError(f"{path}:{line_no}: invalid {label} alignment score: {current.score!r}") from exc
    return Site(row[0], row[1], row[2], mirna, target, orientation)


def read_run(path: Path) -> dict[str, Site]:
    result: dict[str, Site] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for line_no, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            site = parse_site(row, path, line_no)
            if site.read_id in result:
                raise ValueError(f"{path}:{line_no}: duplicate read id {site.read_id!r}")
            result[site.read_id] = site
    if not result:
        raise ValueError(f"{path}: no Hyb rows")
    return result


def read_expression(path: Path | None) -> dict[str, tuple[float, str, str]]:
    if path is None:
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"mirna_id", "expression_value", "expression_unit", "expression_source"}
        expected = ["mirna_id", "expression_value", "expression_unit", "expression_source"]
        if reader.fieldnames != expected:
            raise ValueError(f"{path}: expression header must be exactly {sorted(required)}")
        result = {}
        for line_no, row in enumerate(reader, 2):
            name = row["mirna_id"]
            if not name or name in result:
                raise ValueError(f"{path}:{line_no}: empty or duplicate mirna_id")
            try:
                value = float(row["expression_value"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{path}:{line_no}: invalid expression_value: {row['expression_value']!r}"
                ) from exc
            if not math.isfinite(value):
                raise ValueError(
                    f"{path}:{line_no}: expression_value must be finite: {row['expression_value']!r}"
                )
            result[name] = (value, row["expression_unit"], row["expression_source"])
    return result


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hyb", nargs="+", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--excluded", type=Path, required=True)
    parser.add_argument("--support", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--expression", type=Path)
    parser.add_argument("--expression-threshold", type=float)
    parser.add_argument("--hyb-commit", required=True); parser.add_argument("--hyb-db", required=True)
    parser.add_argument("--run-id", required=True); args = parser.parse_args()
    if len(args.hyb) < 2: parser.error("at least two Hyb runs are required")
    if (args.expression is None) != (args.expression_threshold is None):
        parser.error("expression and expression-threshold must be supplied together")
    if args.expression_threshold is not None and not math.isfinite(args.expression_threshold):
        parser.error("expression threshold must be finite")

    runs = [read_run(path) for path in args.hyb]
    expression = read_expression(args.expression)
    exclusions: list[list[object]] = []; accepted: list[tuple[Site, tuple[float, str, str] | None]] = []
    support_rows: list[list[object]] = []
    all_ids = sorted(set().union(*(run.keys() for run in runs)))
    for read_id in all_ids:
        present = [run.get(read_id) for run in runs]
        observed = [site for site in present if site is not None]
        assignment_counts = Counter(site.assignment_key() for site in observed)
        assignment_examples = {site.assignment_key(): site for site in observed}
        accepted_key = next(iter(assignment_counts)) if len(observed) == len(runs) and len(assignment_counts) == 1 else None
        for rank, (key, count) in enumerate(
            sorted(assignment_counts.items(), key=lambda item: (-item[1], item[0])), 1
        ):
            site = assignment_examples[key]
            supporting_runs = ",".join(
                str(index) for index, candidate in enumerate(present, 1)
                if candidate is not None and candidate.assignment_key() == key
            )
            support_rows.append([
                read_id, rank, site.sequence, site.energy, site.mirna.name, site.target.name,
                site.mirna.read_start, site.mirna.read_end, site.mirna.reference_start,
                site.mirna.reference_end, site.mirna.score, site.target.read_start,
                site.target.read_end, site.target.reference_start, site.target.reference_end,
                site.target.score, count, len(observed), len(runs), supporting_runs,
                "true" if key == accepted_key else "false",
            ])
        if any(site is None for site in present):
            exclusions.append([read_id, "missing_from_replicate", sum(site is not None for site in present), len(runs)]); continue
        sites = [site for site in present if site is not None]
        if len({site.assignment_key() for site in sites}) != 1:
            exclusions.append([read_id, "unstable_assignment", len(sites), len(runs)]); continue
        site = sites[0]; expr = expression.get(site.mirna.name) if expression else None
        if expression and expr is None:
            exclusions.append([read_id, "expression_missing", len(sites), len(runs)]); continue
        if expr is not None and expr[0] < args.expression_threshold:
            exclusions.append([read_id, "below_expression_threshold", len(sites), len(runs)]); continue
        accepted.append((site, expr))
    if not accepted: raise SystemExit("consensus policy retained zero miRNA-mRNA assignments")

    site_rows: list[list[object]] = []; aggregate: dict[tuple[str, str], list[Site]] = defaultdict(list)
    for site, expr in sorted(accepted, key=lambda item: (item[0].mirna.name, item[0].target.name, item[0].read_id)):
        expression_fields: tuple[object, object, object] = ("", "", "") if expr is None else expr
        site_rows.append([
            site.read_id, site.sequence, site.energy, site.mirna.name, site.target.name,
            site.mirna.read_start, site.mirna.read_end, site.mirna.reference_start,
            site.mirna.reference_end, site.mirna.score, site.target.read_start,
            site.target.read_end, site.target.reference_start, site.target.reference_end,
            site.target.score, site.orientation, len(runs), *expression_fields,
        ]); aggregate[(site.mirna.name, site.target.name)].append(site)
    target_rows = []
    for (mirna, target), sites in sorted(aggregate.items()):
        unique_sites = len({(s.target.reference_start, s.target.reference_end) for s in sites})
        target_rows.append([mirna, target, len(sites), unique_sites])

    write_tsv(args.sites, SITE_HEADER, site_rows)
    write_tsv(args.targets, ["mirna_id", "target_id", "consensus_read_ids", "unique_target_sites"], target_rows)
    write_tsv(args.excluded, ["read_id", "reason", "runs_present", "runs_required"], exclusions)
    write_tsv(args.support, [
        "read_id", "assignment_rank", "sequence", "folding_energy", "mirna_id", "target_id",
        "mirna_read_start", "mirna_read_end", "mirna_reference_start", "mirna_reference_end",
        "mirna_alignment_score", "target_read_start", "target_read_end", "target_reference_start",
        "target_reference_end", "target_alignment_score", "runs_supporting_assignment", "runs_present",
        "runs_required", "supporting_run_indices", "accepted",
    ], support_rows)
    reason_counts = Counter(row[1] for row in exclusions)
    manifest = {
        "schema_version": "ago-clip-consensus-1", "status": "complete",
        "policy": "retain only identical normalized miRNA-mRNA assignments present in every clean Hyb run",
        "hyb": {"commit": args.hyb_commit, "database": args.hyb_db, "goal": "detect", "type": "mim", "threads": 1},
        "determinism_controls": {
            "locale": "C", "timezone": "UTC", "python_hash_seed": "0",
            "perl_hash_seed": "0", "perl_perturb_keys": "0",
            "omp_threads": 1, "openblas_threads": 1, "mkl_threads": 1,
        },
        "run_id": args.run_id, "replicates": len(runs),
        "source_reads": {"name": args.reads.name, "sha256": sha256(args.reads)},
        "inputs": [
            {"path": f"raw/replicate-{index}/{path.name}", "sha256": sha256(path), "rows": len(run)}
            for index, (path, run) in enumerate(zip(args.hyb, runs), 1)
        ],
        "expression": None if args.expression is None else {"name": args.expression.name, "sha256": sha256(args.expression), "threshold": args.expression_threshold},
        "accepted_rows": len(site_rows), "target_rows": len(target_rows), "excluded_rows": len(exclusions),
        "exclusion_reasons": dict(sorted(reason_counts.items())),
        "outputs": {
            "sites.tsv": {"sha256": sha256(args.sites), "rows": len(site_rows)},
            "targets.tsv": {"sha256": sha256(args.targets), "rows": len(target_rows)},
            "excluded.tsv": {"sha256": sha256(args.excluded), "rows": len(exclusions)},
            "support.tsv": {"sha256": sha256(args.support), "rows": len(support_rows)},
        },
    }
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__": main()
