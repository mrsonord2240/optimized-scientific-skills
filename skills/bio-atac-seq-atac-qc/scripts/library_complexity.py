#!/usr/bin/env python3
"""
Library complexity: NRF, PBC1, PBC2 on fragments (paired-end) or 5' read positions (single-end).
Reference: pysam 0.22+

Usage:
    python library_complexity.py sample.bam [--mode auto|paired|single]
        [--min-mapq 30] [--exclude-contigs chrM,MT,M,chrMT]

Input: duplicate-retained (pre-dedup) BAM. Flagged duplicates are counted, not skipped.
Paired mode keys each proper pair once by (chrom, fragment start, fragment end), taken from
read 1 (read 1 must pass --min-mapq; mate MAPQ is not inspected). Single mode keys reads by
(chrom, 5' position, strand). Prints strict JSON; PBC2 is null when no position has 2 reads.
Exits 1 when no read passes the filters.
"""

import argparse
import json
import sys
from collections import Counter

import pysam


def library_complexity(bam, mode='auto', min_mapq=30, exclude=('chrM', 'MT', 'M', 'chrMT')):
    exclude = set(exclude)
    counts = Counter()
    n = Counter()
    with pysam.AlignmentFile(bam, 'rb') as bf:
        for r in bf.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if mode == 'auto':
                mode = 'paired' if r.is_paired else 'single'
            if r.reference_name in exclude:
                n['excluded_contig'] += 1
                continue
            if r.mapping_quality < min_mapq:
                n['below_mapq'] += 1
                continue
            if mode == 'paired':
                if not (r.is_paired and r.is_read1 and r.is_proper_pair) or r.template_length == 0:
                    continue
                start = min(r.reference_start, r.next_reference_start)
                key = (r.reference_name, start, start + abs(r.template_length))
            else:
                pos = r.reference_end if r.is_reverse else r.reference_start
                key = (r.reference_name, pos, r.is_reverse)
            if r.is_duplicate:
                n['flagged_duplicate'] += 1
            counts[key] += 1

    total = sum(counts.values())
    distinct = len(counts)
    hist = Counter(counts.values())
    n1, n2 = hist.get(1, 0), hist.get(2, 0)
    out = {
        'mode': mode,
        'unit': 'fragments' if mode == 'paired' else 'reads',
        'min_mapq': min_mapq,
        'total': total,
        'distinct': distinct,
        'NRF': distinct / total if total else None,
        'PBC1': n1 / distinct if distinct else None,
        'PBC2': n1 / n2 if n2 else None,
        'excluded_contig_records': n['excluded_contig'],
        'below_mapq_records': n['below_mapq'],
        'flagged_duplicates_counted': n['flagged_duplicate'],
    }
    if total and n['flagged_duplicate'] == 0 and n1 == distinct:
        out['warning'] = 'every position is unique: BAM looks deduplicated; NRF/PBC are trivially 1 and meaningless'
    return out


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('bam')
    p.add_argument('--mode', choices=['auto', 'paired', 'single'], default='auto')
    p.add_argument('--min-mapq', type=int, default=30)
    p.add_argument('--exclude-contigs', default='chrM,MT,M,chrMT',
                   help='comma-separated contigs to drop (default: mitochondrial names); "" keeps all')
    a = p.parse_args()
    res = library_complexity(a.bam, a.mode, a.min_mapq, [c for c in a.exclude_contigs.split(',') if c])
    print(json.dumps(res, indent=2, allow_nan=False))
    if res['total'] == 0:
        sys.exit('no reads passed the filters')
