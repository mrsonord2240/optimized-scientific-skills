#!/usr/bin/env python3
'''Fetch reads from a region of an indexed BAM or CRAM (indexes first if the index is missing or stale)'''
# Checked on pysam 0.24.1, samtools 1.24 | Verify API if version differs

import argparse
import os
import re
import sys

import pysam


def index_candidates(path):
    '''Index files htslib looks for, in the order it tries them.'''
    stem = os.path.splitext(path)[0]
    if path.endswith('.cram'):
        return [path + '.crai', stem + '.crai']
    return [path + '.csi', stem + '.csi', path + '.bai', stem + '.bai']


def ensure_indexed(path):
    '''Index if no index exists or any index is older than the file.

    Stale indices are deleted first: htslib prefers .csi over .bai, so a stale .csi
    left beside a fresh .bai would still be the one used. CSI is kept if one existed.
    '''
    found = [i for i in index_candidates(path) if os.path.exists(i)]
    if found and all(os.path.getmtime(i) >= os.path.getmtime(path) for i in found):
        return
    csi = any(i.endswith('.csi') for i in found)
    for i in found:
        os.remove(i)
    print(f'Indexing {path}...', file=sys.stderr)
    pysam.index(*(['-c'] if csi else []), path)


def parse_region(region, references):
    '''samtools-style region -> (contig, 0-based start or None, end or None).

    Accepts 'chr1', 'chr1:1,000-2,000', 'chr1:1000-', and contig names containing ':'
    ('HLA-A*01:01:01:01:1-1500', or braced '{HLA-A*01:01:01:01}:1-1500').
    pysam's own fetch(region=...) rejects all of these.
    '''
    region = region.replace(',', '')
    braced = re.fullmatch(r'\{(.+)\}(?::(.*))?', region)
    if braced:
        contig, coords = braced.group(1), braced.group(2) or ''
    elif region in references:
        contig, coords = region, ''
    else:
        contig, _, coords = region.rpartition(':')
    if contig not in references:
        raise ValueError(f'unknown contig {contig!r} (names are case-sensitive; chr1 != 1)')
    start, _, end = coords.partition('-')
    try:
        return contig, (int(start) - 1 if start else None), (int(end) if end else None)
    except ValueError:
        raise ValueError(f'malformed coordinates {coords!r}') from None


def fetch_region(path, region, reference=None):
    ensure_indexed(path)
    mode = 'rc' if path.endswith('.cram') else 'rb'

    with pysam.AlignmentFile(path, mode, reference_filename=reference) as aln:
        try:
            contig, start, end = parse_region(region, aln.references)
        except ValueError as e:
            sys.exit(f'Bad region {region!r}: {e}')
        count = 0
        try:
            for read in aln.fetch(contig, start, end):
                count += 1
                strand = '-' if read.is_reverse else '+'
                print(f'{read.query_name}\t{read.reference_start + 1}\t{strand}')
        except OSError as e:
            hint = ' -- a CRAM needs its reference: pass --reference ref.fa' if mode == 'rc' and not reference else ''
            sys.exit(f'Read failed: {e}{hint}')

    print(f'\nTotal reads in {region}: {count}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('alignment', help='coordinate-sorted BAM or CRAM')
    ap.add_argument('region', help='e.g. chr1:1000000-2000000 (1-based, inclusive)')
    ap.add_argument('--reference', help='reference FASTA (required for CRAM unless its @SQ UR/REF_PATH resolves)')
    args = ap.parse_args()
    fetch_region(args.alignment, args.region, args.reference)
