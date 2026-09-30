#!/usr/bin/env python3
"""Estimate the mono-nucleosome fragment mode (approximate NRL) from ATAC fragment sizes.

Usage: estimate_nrl.py sample.bam

Collects proper-pair read1 template lengths (0-1500 bp), histograms them in 5 bp
bins, finds density peaks (at least 50 bp apart) and reports the tallest peak in
the 150-250 bp window. Exits with a message when no such peak exists (no
nucleosomal periodicity). This is the mono-fragment mode, which tracks the NRL
only approximately; for precision use autocorrelation on cumulative cleavage coverage.
"""
import sys

import numpy as np
import pysam
from scipy.signal import find_peaks

BIN = 5
MIN_SEPARATION_BP = 50


def estimate_nrl(bam_path):
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        frag_lengths = [abs(r.template_length) for r in bam.fetch()
                        if r.is_proper_pair and r.is_read1 and 0 < abs(r.template_length) < 1500]
    if not frag_lengths:
        raise ValueError('no proper-pair read1 fragments found; is the BAM paired-end and indexed?')
    hist, edges = np.histogram(frag_lengths, bins=1500 // BIN, range=(0, 1500))
    # distance is in histogram bins, not bp
    peaks, _ = find_peaks(hist, distance=MIN_SEPARATION_BP // BIN, prominence=hist.max() * 0.05)
    centres = edges[peaks] + BIN / 2
    in_window = (centres > 150) & (centres < 250)
    if not in_window.any():
        raise ValueError('no fragment-size peak at 150-250 bp (peaks found: %s); '
                         'no nucleosomal periodicity detected' % centres.round().tolist())
    return float(centres[in_window][np.argmax(hist[peaks][in_window])])


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    try:
        print(f'Mono-nucleosome fragment mode (approximate NRL): {estimate_nrl(sys.argv[1]):.0f} bp')
    except ValueError as e:
        sys.exit(f'estimate_nrl: {e}')
