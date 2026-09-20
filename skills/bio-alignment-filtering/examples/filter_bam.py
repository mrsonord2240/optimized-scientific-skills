#!/usr/bin/env python3
'''Filter BAM file by various criteria'''
# Reference: pysam 0.22+, samtools 1.19+ | Verify API if version differs
# Checked on pysam 0.24.1 / samtools 1.24.

import argparse
import re
import sys

import pysam


def parse_region(region, references):
    '''samtools-style region (1-based, inclusive) -> (contig, 0-based start, end) for fetch().

    Accepts contig, contig:start, contig:start-end and commas in numbers. A contig name that
    itself contains a colon is found by matching the whole string against the header first.
    '''
    if region in references:
        return region, None, None
    region = region.replace(',', '')
    if region in references:
        return region, None, None
    match = re.match(r'^(.+):(\d+)(?:-(\d*))?$', region)
    if not match or match.group(1) not in references:
        raise ValueError(f'region {region!r} is not contig, contig:start or contig:start-end '
                         f'with a contig from the BAM header')
    start = max(int(match.group(2)) - 1, 0)
    end = int(match.group(3)) if match.group(3) else None
    return match.group(1), start, end


def filter_bam(input_bam, output_bam, min_mapq=0, remove_duplicates=False,
               primary_only=False, proper_pair=False, region=None):

    with pysam.AlignmentFile(input_bam, 'rb') as infile:
        coordinate_sorted = infile.header.to_dict().get('HD', {}).get('SO') == 'coordinate'
        if region:
            if not infile.has_index():
                sys.exit(f'{input_bam} has no index; run samtools index first (needed for --region)')
            contig, start, end = parse_region(region, infile.references)
            iterator = infile.fetch(contig, start, end)
        else:
            iterator = infile

        kept = 0
        removed = 0
        seen_duplicate_flag = False

        with pysam.AlignmentFile(output_bam, 'wb', header=infile.header) as outfile:
            for read in iterator:
                seen_duplicate_flag |= read.is_duplicate
                if read.is_unmapped:
                    removed += 1
                    continue
                if read.mapping_quality < min_mapq:
                    removed += 1
                    continue
                if remove_duplicates and read.is_duplicate:
                    removed += 1
                    continue
                if primary_only and (read.is_secondary or read.is_supplementary):
                    removed += 1
                    continue
                if proper_pair and not read.is_proper_pair:
                    removed += 1
                    continue

                outfile.write(read)
                kept += 1

    print(f'Kept: {kept:,}')
    print(f'Removed: {removed:,}')
    if remove_duplicates and not seen_duplicate_flag:
        print('WARNING: no read carries the duplicate flag (0x400), so --remove-duplicates '
              'removed nothing. Mark duplicates first (see duplicate-handling).', file=sys.stderr)
    return coordinate_sorted


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter BAM file')
    parser.add_argument('input', help='Input BAM file')
    parser.add_argument('output', help='Output BAM file')
    parser.add_argument('-q', '--min-mapq', type=int, default=0, help='Minimum MAPQ')
    parser.add_argument('-d', '--remove-duplicates', action='store_true',
                        help='Remove reads flagged as duplicates (they must already be marked)')
    parser.add_argument('-p', '--primary-only', action='store_true', help='Primary alignments only')
    parser.add_argument('-P', '--proper-pair', action='store_true', help='Properly paired only')
    parser.add_argument('-r', '--region',
                        help='Region, samtools syntax (1-based, inclusive): chr, chr:start, chr:start-end')

    args = parser.parse_args()

    try:
        coordinate_sorted = filter_bam(args.input, args.output, args.min_mapq, args.remove_duplicates,
                                       args.primary_only, args.proper_pair, args.region)
    except ValueError as err:
        sys.exit(f'error: {err}')

    if coordinate_sorted:
        print('\nIndexing output...')
        pysam.index(args.output)
    else:
        print('\nInput header is not SO:coordinate; output left unindexed (sort it first).')
