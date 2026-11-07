#!/usr/bin/env python3
"""Build a strand-aware fragment-size-by-position V-plot around feature positions.

Usage: vplot.py sample.bam features.bed [vplot.png]

Each BED row's start column is the feature center (TSS or motif center). If a
strand column (6) is present, minus-strand features are mirrored so that
downstream is always to the right. Proper-pair read1 fragments within +/-1000 bp
are counted once each at their true center (leftmost mate start + size/2) into a
(fragment size x position) grid, divided by the number of features used. The
image is drawn on a log1p scale. Rows on chromosomes absent from the BAM are skipped.
"""
import sys

import numpy as np
import pysam


def vplot(bam_path, regions_bed, max_size=600, flank=1000):
    """Return (grid, n_features); grid[size, x] with x = position + flank."""
    grid = np.zeros((max_size, 2 * flank))
    n_features = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        known = set(bam.references)
        with open(regions_bed) as fh:
            for line in fh:
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 2 or cols[0].startswith(('#', 'track', 'browser')):
                    continue
                chrom, center = cols[0], int(cols[1])
                if chrom not in known:
                    continue
                minus = len(cols) > 5 and cols[5] == '-'
                n_features += 1
                for r in bam.fetch(chrom, max(0, center - flank), center + flank):
                    if not r.is_proper_pair or not r.is_read1:
                        continue
                    size = abs(r.template_length)
                    if size <= 0 or size >= max_size:
                        continue
                    left = min(r.reference_start, r.next_reference_start)
                    offset = left + size // 2 - center
                    x = (-offset if minus else offset) + flank
                    if 0 <= x < 2 * flank:
                        grid[size, x] += 1
    if n_features == 0:
        raise ValueError('no BED rows on chromosomes present in the BAM header')
    return grid / n_features, n_features


if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    out = sys.argv[3] if len(sys.argv) == 4 else 'vplot.png'
    try:
        g, n = vplot(sys.argv[1], sys.argv[2])
    except ValueError as e:
        sys.exit(f'vplot: {e}')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.imshow(np.log1p(g * 1000), aspect='auto', origin='lower', cmap='magma',
               extent=[-1000, 1000, 0, 600])
    plt.colorbar(label='log(1 + fragments per 1000 features)')
    plt.xlabel('Distance from feature (bp, strand-oriented)')
    plt.ylabel('Fragment size (bp)')
    plt.title(f'V-plot, {n} features')
    plt.savefig(out, dpi=200, bbox_inches='tight')
