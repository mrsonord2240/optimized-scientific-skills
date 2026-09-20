#!/usr/bin/env python3
"""Count reads that still start (5') or end (3') inside a primer footprint.

Independent of samtools ampliconclip: reads the BAM with pysam and the primer BED directly.
A read counts as residual when its aligned 5' end lies inside a primer of the matching strand
(forward read + '+' primer, reverse read + '-' primer), or when its aligned 3' end lies inside
a primer of the opposite strand. Aligned = M/=/X bases, so a soft-clipped primer no longer counts.

Usage: check_primer_residual.py clipped.bam primers.bed [--three-prime]
Both counts are always printed. Exit 1 when any 5' end remains (with --three-prime, also any
3' end); exit 0 otherwise, 2 on bad input.
Checked on pysam 0.24.1, samtools 1.24.
"""
import sys

import pysam


def load_primers(path):
    plus, minus = [], []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[5] not in ("+", "-"):
                sys.exit(f"BED needs >= 6 tab-separated columns with strand in column 6; got: {line.strip()[:80]}")
            (plus if f[5] == "+" else minus).append((f[0], int(f[1]), int(f[2])))
    return plus, minus


def inside(intervals, contig, pos):
    return any(c == contig and s <= pos < e for c, s, e in intervals)


def main(argv):
    args = [a for a in argv if not a.startswith("--")]
    three = "--three-prime" in argv
    if len(args) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    bam, bed = args
    plus, minus = load_primers(bed)
    n = five = three_n = 0
    with pysam.AlignmentFile(bam) as fh:
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.reference_end is None:
                continue
            n += 1
            c, first, last = r.reference_name, r.reference_start, r.reference_end - 1
            if r.is_reverse:
                five += inside(minus, c, last)
                three_n += inside(plus, c, first)
            else:
                five += inside(plus, c, first)
                three_n += inside(minus, c, last)
    if n == 0:
        print("no mapped reads", file=sys.stderr)
        return 2
    print(f"mapped reads={n}  5' end inside a primer={five} ({100 * five / n:.1f}%)  "
          f"3' end inside a primer={three_n} ({100 * three_n / n:.1f}%)"
          + ("" if three else "  (3' not enforced)"))
    return 1 if five or (three and three_n) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
