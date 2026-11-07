#!/usr/bin/env python3
"""Duplicate rate of a duplicate-marked BAM (primary alignments only, the same denominator as
`samtools view -c -f 1024 -F 2304` over `samtools view -c -F 2304`).

Usage:   python dup_rate.py marked.bam
Needs:   pysam (Linux/macOS or WSL). Checked on pysam 0.24.1.
"""
import sys

import pysam

if len(sys.argv) != 2:
    sys.exit(__doc__)

with pysam.AlignmentFile(sys.argv[1], 'rb') as bam:
    total = 0
    duplicates = 0
    for read in bam:
        if read.is_secondary or read.is_supplementary:
            continue
        total += 1
        if read.is_duplicate:
            duplicates += 1

    print(f'Total: {total}')
    print(f'Duplicates: {duplicates}')
    print(f'Rate: {duplicates/total*100:.2f}%' if total else 'Rate: n/a (no primary alignments)')
