#!/usr/bin/env python3
'''Yield each read overlapping any of several regions exactly once (pysam).

Inputs: an indexed BAM/CRAM and a BED file (0-based, half-open) of regions.
Usage:  python fetch_regions.py input.bam regions.bed [reference.fa]   # prints one SAM line per read
        from fetch_regions import fetch_regions                          # in Python, regions = [(contig, start, end), ...]
Checked on pysam 0.24.1, samtools 1.24: equals `samtools view -M -L regions.bed input.bam`.
'''
import sys

import pysam


def fetch_regions(bam, regions):
    '''Yield each read overlapping any (contig, start, end) region once; 0-based, half-open.'''
    merged = {}
    for contig, start, end in sorted(regions):
        ivs = merged.setdefault(contig, [])
        if ivs and start <= ivs[-1][1]:
            ivs[-1][1] = max(ivs[-1][1], end)
        else:
            ivs.append([start, end])
    for contig, ivs in merged.items():
        prev_end = None
        for start, end in ivs:
            for read in bam.fetch(contig, start, end):
                if prev_end is None or read.reference_start >= prev_end:  # else already yielded for the previous interval
                    yield read
            prev_end = end


def read_bed(path):
    regions = []
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith(('#', 'track', 'browser')):
                f = line.split('\t')
                regions.append((f[0], int(f[1]), int(f[2])))
    return regions


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: fetch_regions.py <indexed.bam|cram> <regions.bed> [reference.fa]', file=sys.stderr)
        sys.exit(1)
    reference = sys.argv[3] if len(sys.argv) > 3 else None
    try:
        with pysam.AlignmentFile(sys.argv[1], reference_filename=reference) as bam:
            for read in fetch_regions(bam, read_bed(sys.argv[2])):
                sys.stdout.write(read.to_string() + '\n')
    except (OSError, ValueError) as e:
        print(f'Error: {e}', file=sys.stderr)
        sys.exit(1)
