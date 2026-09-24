#!/usr/bin/env python3
"""Create a temporary BAM and check the shipped pysam helpers."""

import os
import tempfile

import pysam

from pileup_helpers import allele_counts, allele_frequency, find_variants, pileup_text


def write_read(out, name, sequence):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = sequence
    read.flag = 0
    read.reference_id = 0
    read.reference_start = 0
    read.mapping_quality = 60
    read.cigarstring = f"{len(sequence)}M"
    read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    out.write(read)


with tempfile.TemporaryDirectory() as directory:
    fasta = os.path.join(directory, "reference.fa")
    bam = os.path.join(directory, "reads.bam")
    with open(fasta, "w", newline="\n") as handle:
        handle.write(">chrTest\nACGTACGT\n")
    pysam.faidx(fasta)
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrTest", "LN": 8}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as output:
        write_read(output, "ref", "ACGTACGT")
        write_read(output, "alt", "ACTTACGT")
        write_read(output, "unknown", "ACNTACGT")
    pysam.index(bam)

    assert allele_counts(bam, "chrTest", 2) == {"G": 1, "T": 1}
    assert allele_frequency(bam, "chrTest", 2) == {"G": 0.5, "T": 0.5}
    calls = find_variants(bam, fasta, "chrTest", 0, 8, min_depth=2, min_alt_freq=0.4)
    assert calls == [{"chrom": "chrTest", "pos": 3, "ref": "G", "alt": "T",
                      "depth": 2, "alt_count": 1, "freq": 0.5}]
    rows = list(pileup_text(bam, fasta, "chrTest", 0, 8, compute_baq=False))
    assert len(rows) == 8 and rows[2].split("\t")[3] == "3"

print("pileup helper self-test passed")
