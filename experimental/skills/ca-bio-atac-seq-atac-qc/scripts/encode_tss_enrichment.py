#!/usr/bin/env python3
"""
ENCODE-style TSS enrichment: mean signal in the 100 bp window at the TSS / mean signal in the
outer 100 bp at each end of the +/- flank window, on the average TSS profile.
Reference: pyBigWig 0.3+, numpy 1.26+

Usage:
    python encode_tss_enrichment.py <bw_file> <tss_bed> [--flank 2000] [--output score.json]

The score depends on how the bigWig was built (bin size, read representation, normalisation);
see SKILL.md for the recipe. TSS BED: BED6 (strand in column 6). A 1 bp row is used as given;
for a longer (gene) interval the TSS is start for '+' and end-1 for '-'.
Prints JSON with the score and used/skipped TSS counts; exits 1 if no TSS could be scored.
"""

import argparse
import json
import sys

import numpy as np
import pyBigWig


def encode_tss_enrichment(bw_path, tss_bed, flank=2000):
    bw = pyBigWig.open(bw_path)
    chroms = bw.chroms()
    profiles = []
    skipped = {'chrom_not_in_bigwig': 0, 'window_outside_chrom': 0, 'no_data': 0}
    with open(tss_bed, encoding='utf-8') as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6 or parts[5] not in ('+', '-'):
                raise ValueError(f'{tss_bed}:{n}: BED6 with strand (+/-) in column 6 required')
            chrom, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[5]
            tss = end - 1 if strand == '-' else start
            if chrom not in chroms:
                skipped['chrom_not_in_bigwig'] += 1
                continue
            if tss - flank < 0 or tss + flank > chroms[chrom]:
                skipped['window_outside_chrom'] += 1
                continue
            vals = bw.values(chrom, tss - flank, tss + flank)
            if vals is None or len(vals) != 2 * flank:
                skipped['no_data'] += 1
                continue
            vals = np.nan_to_num(np.asarray(vals, dtype=float))
            profiles.append(vals[::-1] if strand == '-' else vals)
    bw.close()

    result = {'TSS_enrichment': None, 'flank': flank, 'tss_used': len(profiles),
              'tss_skipped': skipped}
    if not profiles:
        return result
    avg = np.mean(profiles, axis=0)
    flank_signal = np.mean(np.concatenate([avg[:100], avg[-100:]]))
    center_signal = np.mean(avg[flank - 50: flank + 50])
    if flank_signal > 0:
        result['TSS_enrichment'] = float(center_signal / flank_signal)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Compute ENCODE-style TSS enrichment')
    p.add_argument('bw_file', help='bigWig signal file')
    p.add_argument('tss_bed', help='TSS BED6 (strand in column 6)')
    p.add_argument('--flank', type=int, default=2000, help='flank size (default 2000)')
    p.add_argument('--output', help='output JSON file (default: stdout)')
    a = p.parse_args()
    try:
        res = encode_tss_enrichment(a.bw_file, a.tss_bed, flank=a.flank)
    except ValueError as e:
        sys.exit(f'error: {e}')
    text = json.dumps(res, indent=2, allow_nan=False)
    if a.output:
        with open(a.output, 'w', encoding='utf-8') as f:
            f.write(text)
    print(text)
    if res['TSS_enrichment'] is None:
        sys.exit('error: no TSS scored (check chromosome names, genome build, bigWig coverage)')
