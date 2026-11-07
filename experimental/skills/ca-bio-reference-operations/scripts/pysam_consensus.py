#!/usr/bin/env python3
"""Majority-vote consensus over a BAM window, and its per-position comparison to the reference.

Pedagogical, NOT production: use `samtools consensus` for real work (see references/python-consensus.md).

Inputs:  a sorted, indexed BAM; for `compare`, the FASTA the BAM was aligned to.
         START/END are 0-based half-open (pysam convention); output positions from `compare` are 1-based.
Usage:   pysam_consensus.py consensus BAM CHROM START END [--min-depth 3]      # prints the consensus string
         pysam_consensus.py compare   BAM REF CHROM START END [--min-depth 3]  # prints pos<TAB>ref<TAB>consensus
Tested:  pysam 0.24.1, samtools 1.24
"""
import argparse
from collections import Counter

import pysam


def build_consensus(bam_path, chrom, start, end, min_depth=3):
    """Majority vote over [start, end), 0-based half-open; N where depth < min_depth."""
    consensus = ['N'] * (end - start)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup in bam.pileup(chrom, start, end, truncate=True, max_depth=1_000_000):
            bases = Counter()
            for read in pileup.pileups:
                if not read.is_del and not read.is_refskip:
                    base = read.alignment.query_sequence[read.query_position]
                    bases[base.upper()] += 1

            if sum(bases.values()) >= min_depth:
                consensus[pileup.reference_pos - start] = bases.most_common(1)[0][0]

    return ''.join(consensus)


def compare_to_ref(bam_path, ref_path, chrom, start, end, min_depth=3):
    """[(1-based position, ref base, consensus base)] for called bases that differ from the reference."""
    consensus = build_consensus(bam_path, chrom, start, end, min_depth)
    with pysam.FastaFile(ref_path) as ref:
        reference = ref.fetch(chrom, start, end).upper()   # soft-masked FASTA is lowercase
    return [(start + i + 1, r, c)
            for i, (c, r) in enumerate(zip(consensus, reference))
            if c != 'N' and c != r]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('consensus')
    c.add_argument('bam'); c.add_argument('chrom'); c.add_argument('start', type=int); c.add_argument('end', type=int)
    d = sub.add_parser('compare')
    d.add_argument('bam'); d.add_argument('ref'); d.add_argument('chrom')
    d.add_argument('start', type=int); d.add_argument('end', type=int)
    for p in (c, d):
        p.add_argument('--min-depth', type=int, default=3)
    a = ap.parse_args()
    if a.cmd == 'consensus':
        print(build_consensus(a.bam, a.chrom, a.start, a.end, a.min_depth))
    else:
        for pos, r, cb in compare_to_ref(a.bam, a.ref, a.chrom, a.start, a.end, a.min_depth):
            print(f'{pos}\t{r}\t{cb}')


if __name__ == '__main__':
    main()
