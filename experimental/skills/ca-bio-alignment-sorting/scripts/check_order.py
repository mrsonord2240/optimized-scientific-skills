#!/usr/bin/env python3
"""Check that a BAM/CRAM is really coordinate-sorted by walking its records.

Usage: check_order.py <input.bam|cram> [reference.fa]
Prints the @HD line, then PASS or FAIL. Exit 0 if the records are in (reference order, POS)
order with unmapped reads last, 1 otherwise, 2 for a file with no reference sequences.
The header alone is a hint only: a shuffled BAM labelled SO:coordinate passes a header check.
For CRAM pass the reference (or set REF_PATH).
Checked on pysam 0.24.1.
"""
import sys

import pysam


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 1
    path = sys.argv[1]
    ref = sys.argv[2] if len(sys.argv) > 2 else None

    with pysam.AlignmentFile(path, "rb", check_sq=False, reference_filename=ref) as bam:
        hd = bam.header.to_dict().get("HD", {})
        print("@HD:", " ".join(f"{k}:{v}" for k, v in hd.items()) or "missing")
        if not bam.references:
            print("FAIL: no @SQ lines (uBAM); there is no coordinate order to verify")
            return 2
        prev = (-1, -1)
        n = 0
        for r in bam.fetch(until_eof=True):
            n += 1
            if r.reference_id >= 0:
                key = (r.reference_id, r.reference_start)
            else:
                key = (float("inf"), float("inf"))
            if key < prev:
                print(f"FAIL: record {n} ({r.query_name}) is out of coordinate order")
                return 1
            prev = key
    print(f"PASS: {n} records in coordinate order")
    return 0


if __name__ == "__main__":
    sys.exit(main())
