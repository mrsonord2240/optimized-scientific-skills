#!/usr/bin/env python3
"""Extract a protocol-declared R2 UMI for the pinned targeted miR-eCLIP route."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import tempfile
from itertools import zip_longest
from pathlib import Path
from typing import Iterator, TextIO

UPSTREAM_COMMIT = "75fe74e90e6e4ca670a5af76836d80db09bdbcb1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_fastq(path: Path) -> TextIO:
    with path.open("rb") as handle:
        if handle.read(2) != b"\x1f\x8b":
            raise ValueError(f"{path}: targeted route requires gzip-compressed FASTQ")
    return gzip.open(path, "rt", encoding="utf-8", newline="")


def records(path: Path) -> Iterator[tuple[str, str, str, str]]:
    with open_fastq(path) as handle:
        line_no = 0
        while True:
            header = handle.readline(); line_no += 1
            if not header:
                return
            sequence = handle.readline(); plus = handle.readline(); quality = handle.readline(); line_no += 3
            if not sequence or not plus or not quality:
                raise ValueError(f"{path}:{line_no - 3}: truncated FASTQ record")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"{path}:{line_no - 3}: invalid FASTQ header or separator")
            sequence_value = sequence.rstrip("\r\n")
            quality_value = quality.rstrip("\r\n")
            if len(sequence_value) != len(quality_value):
                raise ValueError(f"{path}:{line_no - 3}: sequence and quality lengths differ")
            yield header, sequence_value, plus, quality_value


def read_id(header: str) -> str:
    return header[1:].strip().split(maxsplit=1)[0]


def output_header(header: str, umi: str) -> str:
    body = header[1:].rstrip("\r\n")
    parts = body.split(maxsplit=1)
    suffix = "" if len(parts) == 1 else " " + parts[1]
    return f"@{parts[0]}_{umi}{suffix}\n"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract an explicitly declared R2-prefix UMI and append it to paired R1 names."
    )
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--output-fastq", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--umi-length", type=int, required=True)
    parser.add_argument("--library-id", required=True)
    parser.add_argument(
        "--protocol-source", required=True,
        help="Versioned protocol, methods section, or library record that declares the UMI length.",
    )
    args = parser.parse_args()

    if not 1 <= args.umi_length <= 64:
        parser.error("--umi-length must be between 1 and 64")
    if not args.library_id.strip() or not args.protocol_source.strip():
        parser.error("--library-id and --protocol-source must be non-empty")
    for path in (args.read1, args.read2):
        if not path.is_file() or path.stat().st_size == 0:
            parser.error(f"missing or empty input: {path}")
    if args.output_fastq.resolve() == args.manifest.resolve():
        parser.error("--output-fastq and --manifest must be different paths")
    for path in (args.output_fastq, args.manifest):
        if path.exists():
            parser.error(f"refusing to overwrite existing output: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)

    output_temp: Path | None = None
    manifest_temp: Path | None = None
    try:
        output_fd, output_name = tempfile.mkstemp(prefix=f".{args.output_fastq.name}.", dir=args.output_fastq.parent)
        os.close(output_fd); output_temp = Path(output_name)
        read_count = 0
        with output_temp.open("w", encoding="utf-8", newline="") as output:
            for index, pair in enumerate(zip_longest(records(args.read1), records(args.read2)), 1):
                r1, r2 = pair
                if r1 is None or r2 is None:
                    raise ValueError(f"paired FASTQ record counts differ at record {index}")
                if read_id(r1[0]) != read_id(r2[0]):
                    raise ValueError(f"paired FASTQ ids differ at record {index}: {read_id(r1[0])!r} != {read_id(r2[0])!r}")
                if len(r2[1]) < args.umi_length:
                    raise ValueError(
                        f"R2 record {index} is shorter than declared UMI length {args.umi_length}"
                    )
                umi = r2[1][:args.umi_length]
                output.write(output_header(r1[0], umi))
                output.write(r1[1] + "\n"); output.write(r1[2]); output.write(r1[3] + "\n")
                read_count += 1
        if read_count == 0:
            raise ValueError("input FASTQs contain no records")

        manifest = {
            "schema_version": "targeted-mir-eclip-umi-1",
            "status": "complete",
            "library_id": args.library_id,
            "protocol_source": args.protocol_source,
            "declared_umi_length": args.umi_length,
            "umi_location": "read2_5prime_prefix",
            "upstream_source_commit": UPSTREAM_COMMIT,
            "upstream_interface_note": (
                "The pinned script's explicit --umi_length is not type-converted; this local route "
                "requires and validates an integer rather than using the conflicting upstream default."
            ),
            "inputs": {
                "read1": {"name": args.read1.name, "sha256": sha256(args.read1)},
                "read2": {"name": args.read2.name, "sha256": sha256(args.read2)},
            },
            "output": {"name": args.output_fastq.name, "sha256": sha256(output_temp), "reads": read_count},
        }
        manifest_fd, manifest_name = tempfile.mkstemp(prefix=f".{args.manifest.name}.", dir=args.manifest.parent)
        os.close(manifest_fd); manifest_temp = Path(manifest_name)
        manifest_temp.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(output_temp, args.output_fastq); output_temp = None
        os.replace(manifest_temp, args.manifest); manifest_temp = None
    finally:
        for path in (output_temp, manifest_temp):
            if path is not None:
                path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
