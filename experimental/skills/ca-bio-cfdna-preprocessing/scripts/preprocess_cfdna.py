#!/usr/bin/env python3
"""cfDNA preprocessing with UMI-aware fgbio consensus.

The executable chain is:

extract UMI -> query-group uBAM -> FASTQ -> align/zip metadata ->
template-coordinate sort -> group -> call consensus -> re-align/zip metadata ->
queryname sort -> filter -> coordinate sort/index.
"""

# Reference: bwa 0.7.17+, fgbio 2.1+, numpy 1.26+, pysam 0.22+, samtools 1.19+
# Verify live interfaces when versions differ.

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pysam


_SAM_TAG = re.compile(r"^[A-Za-z][A-Za-z0-9]$")
_READ_STRUCTURE_SEGMENT = re.compile(r"(?:\d+|\+)([A-Z])")


def _run(cmd: Sequence[str]) -> None:
    """Run one argv-only command and fail on a non-zero exit."""
    subprocess.run([str(part) for part in cmd], check=True)


def _run_pipeline(stages: Iterable[Sequence[str]]) -> None:
    """Connect argv-only processes without invoking a command shell."""
    commands = [[str(part) for part in stage] for stage in stages]
    if len(commands) < 2:
        raise ValueError("a pipeline requires at least two command stages")

    processes: list[subprocess.Popen[bytes]] = []
    previous_stdout = None
    try:
        for command in commands:
            process = subprocess.Popen(
                command,
                stdin=previous_stdout,
                stdout=subprocess.PIPE,
            )
            if previous_stdout is not None:
                previous_stdout.close()
            previous_stdout = process.stdout
            processes.append(process)

        # Consume final stdout even when the final program writes scientific
        # output to an explicit file.
        processes[-1].communicate()
        for process in reversed(processes[:-1]):
            process.wait()
    except BaseException:
        for process in processes:
            if process.poll() is None:
                process.terminate()
        for process in processes:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
        raise

    # Prefer the downstream error: upstream SIGPIPE is commonly a consequence
    # of a later stage rejecting its input.
    for command, process in reversed(list(zip(commands, processes))):
        if process.returncode:
            raise subprocess.CalledProcessError(process.returncode, command)


def _positive_threads(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("threads must be an integer")
    if not 1 <= value <= 256:
        raise ValueError("threads must be between 1 and 256")
    return value


def _input_file(value: os.PathLike[str] | str, label: str) -> Path:
    path = Path(value)
    if not path.is_file():
        raise ValueError(f"{label} must be an existing file: {path}")
    if not os.access(path, os.R_OK):
        raise ValueError(f"{label} is not readable: {path}")
    return path


def _output_file(value: os.PathLike[str] | str, inputs: Sequence[Path]) -> Path:
    path = Path(value)
    if path.exists() and not path.is_file():
        raise ValueError(f"output_bam must be a file path: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    if not os.access(path.parent, os.W_OK):
        raise ValueError(f"output directory is not writable: {path.parent}")
    resolved_output = path.resolve(strict=False)
    if any(resolved_output == source.resolve() for source in inputs):
        raise ValueError("output_bam must differ from every input path")
    return path


def _molecular_segments(read_structures: Sequence[str]) -> int:
    count = 0
    for structure in read_structures:
        segments = _READ_STRUCTURE_SEGMENT.findall(structure)
        if not segments or "T" not in segments:
            raise ValueError(f"invalid read structure (missing template segment): {structure!r}")
        count += segments.count("M")
    return count


def _validated_umi_layout(
    read_structure: str | Sequence[str],
    molecular_index_tags: Sequence[str],
) -> tuple[tuple[str, str], tuple[str, ...]]:
    if isinstance(read_structure, str):
        structures = (read_structure, read_structure)
    else:
        structures = tuple(read_structure)
    if len(structures) != 2 or not all(isinstance(item, str) and item for item in structures):
        raise ValueError("read_structure must define exactly two non-empty mate structures")

    tags = tuple(molecular_index_tags)
    if len(tags) != _molecular_segments(structures):
        raise ValueError("one molecular-index tag is required for every M segment")
    if len(set(tags)) != len(tags) or not all(_SAM_TAG.fullmatch(tag) for tag in tags):
        raise ValueError("molecular_index_tags must be unique two-character SAM tags")
    return (structures[0], structures[1]), tags


def preprocess_cfdna(
    input_bam,
    output_bam,
    reference,
    read_structure="6M11S+T",
    duplex=False,
    threads=8,
    molecular_index_tags=("ZA", "ZB"),
):
    """Build simplex or duplex consensus from a paired, UMI-bearing uBAM.

    ``read_structure`` may be one string applied to both mates or a two-item
    sequence. Each ``M`` segment needs a corresponding two-character tag in
    ``molecular_index_tags``. ``RX`` is also retained as the combined UMI.

    Paths are passed only as argv entries, so spaces and shell metacharacters
    are literal filename characters. The reference must have the indexes and
    sequence dictionary required by bwa, samtools, and fgbio.
    """
    threads = _positive_threads(threads)
    input_bam = _input_file(input_bam, "input_bam")
    reference = _input_file(reference, "reference")
    output_bam = _output_file(output_bam, (input_bam, reference))
    read_structures, molecular_index_tags = _validated_umi_layout(
        read_structure, molecular_index_tags
    )

    work = output_bam.parent
    stem = output_bam.stem
    with_umis = work / f"{stem}_umis.bam"
    with_umis_queryname = work / f"{stem}_umis.queryname.bam"
    aligned_unsorted = work / f"{stem}_aligned.unsorted.bam"
    aligned = work / f"{stem}_aligned.template-coordinate.bam"
    grouped = work / f"{stem}_grouped.bam"
    consensus_unmapped = work / f"{stem}_consensus.unmapped.bam"
    consensus_mapped_unsorted = work / f"{stem}_consensus.mapped.unsorted.bam"
    consensus_mapped_queryname = work / f"{stem}_consensus.mapped.queryname.bam"
    filtered_queryname = work / f"{stem}_filtered.queryname.bam"

    _run(
        [
            "fgbio",
            "ExtractUmisFromBam",
            "--input",
            input_bam,
            "--output",
            with_umis,
            "--read-structure",
            *read_structures,
            "--molecular-index-tags",
            *molecular_index_tags,
            "--single-tag",
            "RX",
        ]
    )

    # ZipperBams requires identically query-grouped unmapped and mapped reads.
    _run(
        [
            "samtools",
            "sort",
            "-n",
            "--threads",
            str(threads),
            "-o",
            with_umis_queryname,
            with_umis,
        ]
    )
    _run_pipeline(
        [
            ["samtools", "fastq", "-T", ",".join(("RX", *molecular_index_tags)), with_umis_queryname],
            [
                "bwa",
                "mem",
                "-C",
                "-t",
                str(threads),
                "-K",
                "150000000",
                "-Y",
                "-p",
                reference,
                "-",
            ],
            [
                "fgbio",
                "ZipperBams",
                "--unmapped",
                with_umis_queryname,
                "--ref",
                reference,
                "--output",
                aligned_unsorted,
            ],
        ]
    )
    _run(
        [
            "samtools",
            "sort",
            "--template-coordinate",
            "--threads",
            str(threads),
            "-o",
            aligned,
            aligned_unsorted,
        ]
    )

    strategy = "paired" if duplex else "adjacency"
    _run(
        [
            "fgbio",
            "GroupReadsByUmi",
            "--input",
            aligned,
            "--output",
            grouped,
            "--strategy",
            strategy,
            "--edits",
            "1",
        ]
    )

    caller = "CallDuplexConsensusReads" if duplex else "CallMolecularConsensusReads"
    _run(
        [
            "fgbio",
            caller,
            "--input",
            grouped,
            "--output",
            consensus_unmapped,
            "--min-reads",
            "1",
        ]
    )

    _run_pipeline(
        [
            ["samtools", "fastq", consensus_unmapped],
            [
                "bwa",
                "mem",
                "-t",
                str(threads),
                "-K",
                "150000000",
                "-Y",
                "-p",
                reference,
                "-",
            ],
            [
                "fgbio",
                "ZipperBams",
                "--unmapped",
                consensus_unmapped,
                "--ref",
                reference,
                "--tags-to-reverse",
                "Consensus",
                "--tags-to-revcomp",
                "Consensus",
                "--output",
                consensus_mapped_unsorted,
            ],
        ]
    )

    # FilterConsensusReads is template-aware and requires query-grouped input.
    _run(
        [
            "samtools",
            "sort",
            "-n",
            "--threads",
            str(threads),
            "-o",
            consensus_mapped_queryname,
            consensus_mapped_unsorted,
        ]
    )
    min_reads = ["2", "1", "1"] if duplex else ["2"]
    _run(
        [
            "fgbio",
            "FilterConsensusReads",
            "--input",
            consensus_mapped_queryname,
            "--output",
            filtered_queryname,
            "--ref",
            reference,
            "--min-reads",
            *min_reads,
            "--max-read-error-rate",
            "0.025",
            "--max-base-error-rate",
            "0.1",
            "--min-base-quality",
            "40",
            "--reverse-per-base-tags",
        ]
    )
    _run(
        [
            "samtools",
            "sort",
            "--threads",
            str(threads),
            "-o",
            output_bam,
            filtered_queryname,
        ]
    )
    _run(["samtools", "index", "--threads", str(threads), output_bam])
    return output_bam


def insert_size_qc(bam_path, max_size=600):
    """Summarize the primary, PF, nonduplicate proper-pair fragment population.

    The positive-template-length convention counts one observation per pair.
    ``filtered_counts`` reports overlapping flag/reason counts, not a partition.
    Biological interpretation still requires the library chemistry, collection
    conditions, and any physical or in-silico size selection.
    """
    if isinstance(max_size, bool) or not isinstance(max_size, int):
        raise TypeError("max_size must be an integer")
    if max_size <= 0:
        raise ValueError("max_size must be greater than zero")

    reasons = {
        "input_records": 0,
        "unmapped": 0,
        "not_proper_pair": 0,
        "secondary": 0,
        "supplementary": 0,
        "duplicate": 0,
        "qc_fail": 0,
        "nonpositive_template_length": 0,
        "over_max_size": 0,
    }
    sizes: list[int] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            reasons["input_records"] += 1
            reasons["unmapped"] += int(read.is_unmapped)
            reasons["not_proper_pair"] += int(not read.is_proper_pair)
            reasons["secondary"] += int(read.is_secondary)
            reasons["supplementary"] += int(read.is_supplementary)
            reasons["duplicate"] += int(read.is_duplicate)
            reasons["qc_fail"] += int(read.is_qcfail)
            reasons["nonpositive_template_length"] += int(read.template_length <= 0)
            reasons["over_max_size"] += int(read.template_length > max_size)
            if (
                read.is_proper_pair
                and not read.is_unmapped
                and not read.is_secondary
                and not read.is_supplementary
                and not read.is_duplicate
                and not read.is_qcfail
                and 0 < read.template_length <= max_size
            ):
                sizes.append(read.template_length)

    values = np.asarray(sizes, dtype=np.int64)
    result = {
        "n": len(sizes),
        "mode_bp": 0,
        "median_bp": 0.0,
        "short_frac_90_150": 0.0,
        "frac_over_250bp": 0.0,
        "filtered_counts": {**reasons, "accepted": len(sizes)},
        "interpretation_context": [
            "library_chemistry",
            "collection_and_plasma_processing",
            "physical_or_in_silico_size_selection",
        ],
    }
    if values.size:
        result.update(
            {
                "mode_bp": int(np.bincount(values).argmax()),
                "median_bp": float(np.median(values)),
                "short_frac_90_150": float(np.mean((values >= 90) & (values <= 150))),
                "frac_over_250bp": float(np.mean(values > 250)),
            }
        )
    return result


if __name__ == "__main__":
    print("cfDNA preprocessing pipeline")
    print("preprocess_cfdna(in_bam, out_bam, ref, duplex=False) - simplex UMI consensus")
    print("preprocess_cfdna(in_bam, out_bam, ref, duplex=True)  - duplex consensus")
    print("insert_size_qc(bam) - fragment-length QC with explicit filter counts")
