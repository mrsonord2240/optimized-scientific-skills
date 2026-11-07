#!/usr/bin/env python3
"""Mark duplicates with pysam: sort -n, fixmate -m, sort, markdup, index.

Usage:   python pysam_markdup.py input.bam marked.bam
Needs:   pysam (Linux/macOS or WSL; no Windows wheel). Checked on pysam 0.24.1 (samtools 1.24).
Output:  marked.bam plus marked.bam.bai; prints the input/output record counts and the flagged count.
The intermediates go in a temporary directory, and markdup drops no records, so the counts must match.
"""
import os
import sys
import tempfile

import pysam

if len(sys.argv) != 3:
    sys.exit(__doc__)
input_bam, marked_bam = sys.argv[1], sys.argv[2]

with tempfile.TemporaryDirectory() as tmp:
    namesort = os.path.join(tmp, 'namesort.bam')
    fixmate = os.path.join(tmp, 'fixmate.bam')
    coordsort = os.path.join(tmp, 'coordsort.bam')

    # Sort by name
    pysam.sort('-n', '-o', namesort, input_bam)

    # Fixmate
    pysam.fixmate('-m', namesort, fixmate)

    # Sort by coordinate
    pysam.sort('-o', coordsort, fixmate)

    # Mark duplicates
    pysam.markdup(coordsort, marked_bam)

# Index
pysam.index(marked_bam)

n_in = int(pysam.view('-c', input_bam).strip())
n_out = int(pysam.view('-c', marked_bam).strip())
n_dup = int(pysam.view('-c', '-f', '1024', marked_bam).strip())
print(f'records in: {n_in}, out: {n_out}, flagged duplicate: {n_dup}')
if n_in == 0 or n_in != n_out:
    sys.exit('ERROR: record count changed or input empty')
