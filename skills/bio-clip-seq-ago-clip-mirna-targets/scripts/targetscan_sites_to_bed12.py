#!/usr/bin/env python3
"""Map normalized TargetScan UTR-relative sites to strand-safe genomic BED12."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class UtrMap:
    transcript_id: str
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]


def exact_header(reader: csv.DictReader, required: set[str], path: Path) -> None:
    if reader.fieldnames is None or set(reader.fieldnames) != required:
        raise ValueError(f"{path}: header must be exactly {sorted(required)}")


def load_maps(path: Path, args: argparse.Namespace) -> dict[str, UtrMap]:
    required = {"transcript_id", "chrom", "strand", "utr_exon_starts_0", "utr_exon_ends_0", "assembly", "annotation_release", "targetscan_release"}
    result: dict[str, UtrMap] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t"); exact_header(reader, required, path)
        for line_no, row in enumerate(reader, 2):
            if row["assembly"] != args.assembly or row["annotation_release"] != args.annotation_release or row["targetscan_release"] != args.targetscan_release:
                raise ValueError(f"{path}:{line_no}: release metadata does not match command arguments")
            starts = [int(value) for value in row["utr_exon_starts_0"].rstrip(",").split(",")]
            ends = [int(value) for value in row["utr_exon_ends_0"].rstrip(",").split(",")]
            if len(starts) != len(ends) or not starts or row["strand"] not in {"+", "-"}:
                raise ValueError(f"{path}:{line_no}: invalid exon arrays or strand")
            exons = tuple(zip(starts, ends))
            if any(start < 0 or end <= start for start, end in exons) or list(exons) != sorted(exons):
                raise ValueError(f"{path}:{line_no}: exons must be sorted, non-overlapping genomic half-open intervals")
            if any(exons[index][1] > exons[index + 1][0] for index in range(len(exons) - 1)):
                raise ValueError(f"{path}:{line_no}: overlapping UTR exons")
            transcript = row["transcript_id"]
            if not transcript or transcript in result:
                raise ValueError(f"{path}:{line_no}: empty or duplicate transcript_id")
            result[transcript] = UtrMap(transcript, row["chrom"], row["strand"], exons)
    return result


def project(mapping: UtrMap, start_1: int, end_1: int) -> list[tuple[int, int]]:
    if start_1 < 1 or end_1 < start_1:
        raise ValueError("TargetScan coordinates must be 1-based inclusive with end >= start")
    wanted_start, wanted_end = start_1 - 1, end_1
    tx_cursor = 0; blocks: list[tuple[int, int]] = []
    transcript_order = mapping.exons if mapping.strand == "+" else tuple(reversed(mapping.exons))
    for exon_start, exon_end in transcript_order:
        length = exon_end - exon_start; overlap_start = max(wanted_start, tx_cursor); overlap_end = min(wanted_end, tx_cursor + length)
        if overlap_start < overlap_end:
            local_start, local_end = overlap_start - tx_cursor, overlap_end - tx_cursor
            if mapping.strand == "+": block = (exon_start + local_start, exon_start + local_end)
            else: block = (exon_end - local_end, exon_end - local_start)
            blocks.append(block)
        tx_cursor += length
    if sum(end - start for start, end in blocks) != end_1 - start_1 + 1:
        raise ValueError(f"site {start_1}-{end_1} exceeds mapped UTR length {tx_cursor}")
    return sorted(blocks)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sites", type=Path, required=True, help="normalized TargetScan site TSV")
    parser.add_argument("--utr-map", type=Path, required=True, help="versioned spliced UTR map TSV")
    parser.add_argument("--output", type=Path, required=True); parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--assembly", required=True); parser.add_argument("--annotation-release", required=True)
    parser.add_argument("--targetscan-release", required=True); args = parser.parse_args()
    mappings = load_maps(args.utr_map, args)
    required = {"transcript_id", "utr_start_1", "utr_end_1", "mirna_family", "site_type"}
    rows = []; missing = []; converted = 0
    with args.sites.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t"); exact_header(reader, required, args.sites)
        for line_no, row in enumerate(reader, 2):
            mapping = mappings.get(row["transcript_id"])
            if mapping is None:
                missing.append(row["transcript_id"]); continue
            blocks = project(mapping, int(row["utr_start_1"]), int(row["utr_end_1"]))
            chrom_start, chrom_end = blocks[0][0], blocks[-1][1]
            name = "|".join((row["transcript_id"], row["mirna_family"], row["site_type"]))
            sizes = ",".join(str(end - start) for start, end in blocks) + ","
            offsets = ",".join(str(start - chrom_start) for start, _ in blocks) + ","
            rows.append([mapping.chrom, chrom_start, chrom_end, name, 0, mapping.strand, chrom_start, chrom_end, "0", len(blocks), sizes, offsets])
            converted += 1
    if missing:
        raise SystemExit(f"unmapped transcript ids ({len(missing)} rows); first: {missing[0]}")
    if not rows: raise SystemExit("no TargetScan sites converted")
    rows.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
    manifest = {
        "schema_version": "targetscan-utr-to-bed12-1", "status": "complete", "converted_sites": converted,
        "coordinate_contract": "input 1-based inclusive spliced UTR; output 0-based half-open genomic BED12",
        "assembly": args.assembly, "annotation_release": args.annotation_release,
        "targetscan_release": args.targetscan_release, "strand_preserved": True,
    }
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__": main()
