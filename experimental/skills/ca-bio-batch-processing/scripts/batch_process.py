#!/usr/bin/env python3
"""Safe, deterministic batch operations for FASTA and FASTQ files."""

from __future__ import annotations

import csv
import json
import re
import tempfile
from collections import OrderedDict
from itertools import islice
from multiprocessing import get_context
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


SUMMARY_FIELDS = ("file", "sequences", "total_bp", "min_len", "max_len", "avg_len")
_SAFE_PREFIX = re.compile(r"^[A-Za-z0-9][A-Za-z0-9.-]*$")
_WINDOWS_RESERVED = {
    "CON",
    "PRN",
    "AUX",
    "NUL",
    *(f"COM{number}" for number in range(1, 10)),
    *(f"LPT{number}" for number in range(1, 10)),
}


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer; got {value!r}")
    return value


def stable_paths(directory, pattern, recursive=False):
    """Discover files in a reproducible relative-path order."""
    root = Path(directory)
    iterator = root.rglob(pattern) if recursive else root.glob(pattern)

    def stable_key(path):
        relative = path.relative_to(root).as_posix()
        return relative.casefold(), relative

    return sorted((path for path in iterator if path.is_file()), key=stable_key)


def count_streaming(filepath, format):
    """Count records without materializing them."""
    with Path(filepath).open("r", encoding="utf-8") as handle:
        return sum(1 for _ in SeqIO.parse(handle, format))


def split_by_count(input_file, format, records_per_file, output_prefix):
    """Split into positive-size chunks while retaining only one chunk in memory."""
    _positive_integer(records_per_file, "records_per_file")
    written = []
    file_num = 1
    with Path(input_file).open("r", encoding="utf-8") as input_handle:
        records = SeqIO.parse(input_handle, format)
        while True:
            batch = list(islice(records, records_per_file))
            if not batch:
                break
            out = f"{output_prefix}_{file_num}.{format}"
            SeqIO.write(batch, out, format)
            written.append(out)
            file_num += 1
    return written


def _validated_prefix(record_id):
    prefix = record_id.split("_", 1)[0]
    if not _SAFE_PREFIX.fullmatch(prefix):
        raise ValueError(
            f"unsafe record-id prefix {prefix!r}: use an ASCII letter/digit first, "
            "then only letters, digits, dot, or hyphen"
        )
    if prefix.split(".", 1)[0].upper() in _WINDOWS_RESERVED:
        raise ValueError(f"unsafe record-id prefix {prefix!r}: reserved filename")
    return prefix


def _write_manifest_exclusive(path, payload):
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def split_by_prefix(
    input_file,
    format,
    output_directory,
    *,
    max_open_files=64,
    manifest_name="split-prefix-manifest.json",
):
    """Split safely by id prefix with bounded handles and a completion manifest.

    Sequence files and the manifest are created exclusively: an existing target
    is never truncated. If parsing or writing fails, open handles are closed and
    the manifest records ``status: partial`` plus the triggering error.
    """
    _positive_integer(max_open_files, "max_open_files")
    output_root = Path(output_directory).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = (output_root / manifest_name).resolve()
    if manifest_path.parent != output_root:
        raise ValueError("manifest_name must name a file directly inside output_directory")
    if manifest_path.exists():
        raise FileExistsError(f"refusing to overwrite manifest: {manifest_path}")

    handles = OrderedDict()
    targets = {}
    canonical_names = {}
    counts = {}
    records_written = 0
    failure = None
    try:
        with Path(input_file).open("r", encoding="utf-8") as input_handle:
            for record in SeqIO.parse(input_handle, format):
                prefix = _validated_prefix(record.id)
                filename = f"{prefix}.{format}"
                canonical = filename.casefold()
                other = canonical_names.setdefault(canonical, prefix)
                if other != prefix:
                    raise ValueError(
                        f"prefixes {other!r} and {prefix!r} collide on case-insensitive filesystems"
                    )

                handle = handles.pop(prefix, None)
                if handle is None:
                    if len(handles) >= max_open_files:
                        _, oldest = handles.popitem(last=False)
                        oldest.close()
                    if prefix not in targets:
                        target = (output_root / filename).resolve()
                        if target.parent != output_root:
                            raise ValueError(f"derived path escapes output_directory: {target}")
                        handle = target.open("x", encoding="utf-8", newline="\n")
                        targets[prefix] = target
                        counts[prefix] = 0
                    else:
                        handle = targets[prefix].open("a", encoding="utf-8", newline="\n")
                handles[prefix] = handle

                SeqIO.write(record, handle, format)
                counts[prefix] += 1
                records_written += 1
    except Exception as exc:
        failure = {"type": type(exc).__name__, "message": str(exc)}
        raise
    finally:
        for handle in handles.values():
            handle.close()
        payload = {
            "status": "partial" if failure else "complete",
            "input": str(Path(input_file)),
            "format": format,
            "records_written": records_written,
            "files": [
                {
                    "prefix": prefix,
                    "path": str(target.relative_to(output_root)),
                    "records": counts[prefix],
                }
                for prefix, target in sorted(
                    targets.items(), key=lambda item: (item[0].casefold(), item[0])
                )
            ],
            "error": failure,
        }
        _write_manifest_exclusive(manifest_path, payload)
    return payload


def summarize_files(files, format, output_csv):
    """Write deterministic summaries, including a header-only CSV for no files."""
    summaries = []
    for filepath in sorted(
        (Path(path) for path in files),
        key=lambda path: (path.as_posix().casefold(), path.as_posix()),
    ):
        count = total = 0
        min_len = None
        max_len = 0
        with filepath.open("r", encoding="utf-8") as input_handle:
            for record in SeqIO.parse(input_handle, format):
                length = len(record.seq)
                count += 1
                total += length
                max_len = max(max_len, length)
                min_len = length if min_len is None else min(min_len, length)
        summaries.append(
            {
                "file": filepath.name,
                "sequences": count,
                "total_bp": total,
                "min_len": min_len or 0,
                "max_len": max_len,
                "avg_len": total / count if count else 0,
            }
        )

    with Path(output_csv).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(summaries)
    return summaries


def process_file(filepath):
    """Spawn-safe worker: the function is importable at module scope."""
    filepath = Path(filepath)
    total = 0
    bp = 0
    with filepath.open("r", encoding="utf-8") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            total += 1
            bp += len(record.seq)
    return {"file": filepath.name, "count": total, "total_bp": bp}


def process_files_parallel(files, workers=4, start_method=None):
    """Process a stable file list under fork or spawn."""
    _positive_integer(workers, "workers")
    ordered_files = sorted(
        (Path(path) for path in files),
        key=lambda path: (path.as_posix().casefold(), path.as_posix()),
    )
    context = get_context(start_method) if start_method else get_context()
    with context.Pool(workers) as pool:
        return pool.map(process_file, ordered_files)


def make_demo_fastas(directory, n_files, seqs_per_file):
    """Write small demo FASTA files so the example runs standalone."""
    paths = []
    for fi in range(n_files):
        path = Path(directory) / f"sample{fi}.fasta"
        records = [
            SeqRecord(Seq("ACGT" * (10 + i)), id=f"s{fi}_{i}", description="")
            for i in range(seqs_per_file)
        ]
        SeqIO.write(records, path, "fasta")
        paths.append(str(path))
    return paths


def main():
    with tempfile.TemporaryDirectory() as workdir:
        files = make_demo_fastas(workdir, n_files=3, seqs_per_file=5)

        for filepath in stable_paths(workdir, "*.fasta"):
            print(f"{filepath.name}: {count_streaming(filepath, 'fasta')} sequences")

        chunks = split_by_count(files[0], "fasta", 2, str(Path(workdir) / "chunk"))
        print(f"Split {Path(files[0]).name} into {len(chunks)} chunks of <=2 records")

        index_path = str(Path(workdir) / "combined.idx")
        records = SeqIO.index_db(index_path, files, "fasta")
        print(f"Indexed {len(records)} records across {len(files)} files")
        print(f"Random lookup s2_3 length: {len(records['s2_3'].seq)}")
        records.close()


if __name__ == "__main__":
    main()
