#!/usr/bin/env python3
'''Extract the reads overlapping BED regions: each read once, in coordinate order.

Same records as `samtools view -L targets.bed` for tab- or space-delimited rows with
start < end (samtools reads a zero-width row, start == end, by its own rule, which fetch()
cannot express: drop or widen such rows). Needs a coordinate-sorted, indexed BAM.

Usage: python filter_by_bed.py input.bam targets.bed output.bam
Checked on pysam 0.24.1 / samtools 1.24.
'''
import sys

import pysam


def read_bed(bed_path):
    regions = []
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.split()
            regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions


def main(input_bam, bed_path, output_bam):
    with pysam.AlignmentFile(input_bam, 'rb') as infile:
        order = {name: i for i, name in enumerate(infile.references)}
        regions = sorted((r for r in read_bed(bed_path) if r[0] in order),
                         key=lambda r: (order[r[0]], r[1]))
        merged = []
        for chrom, start, end in regions:
            if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
                merged[-1][2] = max(merged[-1][2], end)
            else:
                merged.append([chrom, start, end])

        with pysam.AlignmentFile(output_bam, 'wb', header=infile.header) as outfile:
            prev_chrom, prev_end = None, 0
            for chrom, start, end in merged:
                if chrom != prev_chrom:
                    prev_end = 0
                for read in infile.fetch(chrom, start, end):
                    if read.reference_start >= prev_end:
                        outfile.write(read)
                prev_chrom, prev_end = chrom, end


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
