#!/usr/bin/env python3
'''View SAM/BAM/CRAM contents using pysam'''
# Reference: pysam 0.24.1, samtools 1.24 (checked); pysam 0.22+ | Verify API if version differs

import sys

import pysam


def count_mapped(path, reference):
    '''Return (mapped, unmapped, how). BAM with an index: read the index. SAM, CRAM, unindexed BAM: scan.'''
    try:
        with pysam.AlignmentFile(path, reference_filename=reference) as bam:
            # index counts are only reliable for BAM: pysam returns 0/0 for an indexed CRAM
            if bam.is_bam:
                stats = bam.get_index_statistics()
                mapped = sum(s.mapped for s in stats)
                unmapped = sum(s.unmapped for s in stats) + bam.nocoordinate
                return mapped, unmapped, 'index'
    except (ValueError, AttributeError, OSError):
        pass  # no index: count records by reading the whole file
    mapped = unmapped = 0
    with pysam.AlignmentFile(path, reference_filename=reference) as bam:
        for read in bam:
            if read.is_unmapped:
                unmapped += 1
            else:
                mapped += 1
    return mapped, unmapped, 'scan'


def view_bam(path, limit=10, reference=None):
    # mode defaults to 'r': format (SAM/BAM/CRAM) is detected from the file
    with pysam.AlignmentFile(path, reference_filename=reference) as bam:
        print(f'References: {bam.nreferences}')
    mapped, unmapped, how = count_mapped(path, reference)
    print(f'Mapped: {mapped}  Unmapped: {unmapped}  (from {how})')
    print()
    print('name\tchrom:start(0-based)\tstrand\tcigar')

    with pysam.AlignmentFile(path, reference_filename=reference) as bam:
        for i, read in enumerate(bam):
            if i >= limit:
                break
            strand = '-' if read.is_reverse else '+'
            if read.is_unmapped:
                where = f'{read.reference_name}:{read.reference_start}' if read.reference_name else '*'
                print(f'{read.query_name}\t{where}\t{strand}\tunmapped')
            else:
                print(f'{read.query_name}\t{read.reference_name}:{read.reference_start}\t{strand}\t{read.cigarstring}')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: view_bam.py <input.sam|bam|cram> [limit] [reference.fa]')
        sys.exit(1)

    path = sys.argv[1]
    try:
        limit = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    except ValueError:
        print(f'Error: limit must be an integer, got {sys.argv[2]!r}', file=sys.stderr)
        print('Usage: view_bam.py <input.sam|bam|cram> [limit] [reference.fa]', file=sys.stderr)
        sys.exit(1)
    reference = sys.argv[3] if len(sys.argv) > 3 else None
    try:
        view_bam(path, limit, reference)
    except (OSError, ValueError) as e:
        hint = ' (CRAM needs its reference: pass reference.fa as the 3rd argument)' if path.endswith('.cram') else ''
        print(f'Error reading {path}: {e}{hint}', file=sys.stderr)
        sys.exit(1)
