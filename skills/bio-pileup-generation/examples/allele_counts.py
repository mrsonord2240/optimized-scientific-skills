#!/usr/bin/env python3
'''Count alleles at a single position using pysam pileup'''
# Checked on samtools 1.24, pysam 0.24.1 | Verify API if version differs
# Without a fastafile pysam applies no BAQ, i.e. equals `samtools mpileup -B`; see SKILL.md for the full parameter table.

import sys
import pysam

from pileup_helpers import allele_counts

def main():
    if len(sys.argv) < 3:
        print('Usage: allele_counts.py <input.bam> <contig:position>   (position is 1-based)')
        print('Example: allele_counts.py sample.bam chr1:1000000')
        sys.exit(1)

    bam_path = sys.argv[1]
    region = sys.argv[2].replace(',', '')

    chrom, sep, pos = region.rpartition(':')  # contig names may contain ':' (HLA-A*01:01:01:01)
    if not sep or not pos.isdigit() or int(pos) < 1:
        sys.exit(f'Error: expected <contig>:<1-based position>, got {sys.argv[2]!r} (ranges are not supported)')
    pos = int(pos) - 1  # Convert to 0-based

    try:
        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            if chrom not in bam.references:
                sys.exit(f'Error: contig {chrom!r} is not in the BAM header (see: samtools view -H {bam_path} | grep ^@SQ)')
            if not bam.has_index():
                sys.exit(f'Error: {bam_path} has no index (run: samtools index {bam_path})')
        counts = allele_counts(bam_path, chrom, pos, min_mapping_quality=20,
                               min_base_quality=20)
    except (OSError, ValueError) as e:
        sys.exit(f'Error: cannot read {bam_path}: {e}')
    total = sum(counts.values())

    print(f'Position: {chrom}:{pos+1}')
    print(f'Total depth: {total}')
    print('Allele counts:')

    for base, count in sorted(counts.items(), key=lambda x: -x[1]):
        freq = count / total * 100 if total > 0 else 0
        print(f'  {base}: {count} ({freq:.1f}%)')

if __name__ == '__main__':
    main()
