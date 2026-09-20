# Reference: numpy 1.26+, pysam 0.22+ | Verify API if version differs (checked on pysam 0.24.1, numpy 2.5)
#
# Usage: python validate_alignment.py sample.bam [-n N]
# Exit status: 0 = every metric PASS or WARN, 1 = at least one metric in the FAIL band,
#              2 = file unreadable / truncated / no records (nothing was validated).
# Reads are counted through the whole file with until_eof=True: no index needed, and unplaced
# unmapped reads (RNAME '*', the tail of a coordinate-sorted BAM) are included. Secondary and
# supplementary records are skipped in every rate (same denominators as `samtools flagstat` "primary").
#
# Bands (germline short-read DNA defaults; assay-specific, see SKILL.md "Quality Thresholds Summary"):
#   mapping rate   > 95 PASS, 90-95 WARN, < 90 FAIL  -- ATAC/ChIP/RNA/aDNA can be lower.
#   proper pair    > 90 PASS, 80-90 WARN, < 80 FAIL  -- spliced RNA-seq and amplicon lower; n/a for single-end.
#   strand balance forward fraction 0.48-0.52 PASS, 0.45-0.55 WARN -- autosomal DNA-seq only; deviates for stranded RNA / bisulfite.
#   mean MAPQ      > 40 PASS, 30-40 WARN, < 30 FAIL  -- bimodal, aligner-specific (wrong for STAR sentinel 255).
# -n N samples only the first N records: head-of-file biased, and in a coordinate-sorted BAM the unplaced
# unmapped tail is not reached, so the mapping rate is then an overestimate. Use `samtools flagstat` for the exact rate.
import sys
import argparse
import numpy as np
import pysam


def grade_min(value, good, warn):
    return 'PASS' if value > good else 'WARN' if value >= warn else 'FAIL'


def grade_strand(fraction):
    return 'PASS' if 0.48 <= fraction <= 0.52 else 'WARN' if 0.45 <= fraction <= 0.55 else 'FAIL'


def validate_bam(bam_file, sample_size=0):
    mapped = unmapped = 0
    proper_pair = paired = 0
    forward = reverse = 0
    insert_sizes = []
    mapqs = []
    n_read = 0

    try:
        with pysam.AlignmentFile(bam_file, 'rb', check_sq=False) as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary:
                    continue
                n_read += 1
                if sample_size and n_read > sample_size:
                    n_read -= 1
                    break

                if read.is_unmapped:
                    unmapped += 1
                    continue
                mapped += 1
                if read.is_reverse:
                    reverse += 1
                else:
                    forward += 1
                mapqs.append(read.mapping_quality)

                if read.is_paired:
                    paired += 1
                    if read.is_proper_pair:
                        proper_pair += 1
                        if read.template_length > 0:
                            insert_sizes.append(read.template_length)
    except (OSError, ValueError) as e:
        print(f'ERROR: cannot read {bam_file} completely: {e}', file=sys.stderr)
        return 2

    total = mapped + unmapped
    if total == 0:
        print(f'ERROR: {bam_file} has no primary records; nothing to validate', file=sys.stderr)
        return 2

    scope = f'first {n_read} primary records' if sample_size and n_read == sample_size else f'all {n_read} primary records'
    print(f'=== Alignment Validation ({scope}) ===\n')
    results = []  # (label, grade)

    print('--- Mapping ---')
    map_rate = 100 * mapped / total
    print(f'Mapped: {mapped} ({map_rate:.2f}%)')
    print(f'Unmapped: {unmapped} ({100*unmapped/total:.2f}%)')
    results.append(('Mapping rate', grade_min(map_rate, 95, 90)))

    print('\n--- Pairing ---')
    if paired > 0:
        pair_rate = 100 * proper_pair / paired
        print(f'Properly paired: {proper_pair} ({pair_rate:.2f}% of mapped paired reads)')
        results.append(('Proper pairing', grade_min(pair_rate, 90, 80)))
    else:
        print('No mapped paired reads (single-end or unpaired): proper pairing not graded')

    print('\n--- Insert Size ---')
    if insert_sizes:
        print(f'Median: {np.median(insert_sizes):.0f} bp')
        print(f'Mean: {np.mean(insert_sizes):.0f} bp')
        print(f'Std: {np.std(insert_sizes):.0f} bp')
    else:
        print('No proper pairs with a positive template length')

    if mapped == 0:
        print('\nNo mapped reads: strand balance and MAPQ not computed')
    else:
        print('\n--- Strand Balance ---')
        strand_ratio = forward / (forward + reverse)
        print(f'Forward: {forward}, Reverse: {reverse}')
        print(f'Forward fraction F/(F+R): {strand_ratio:.3f}')
        results.append(('Strand balance', grade_strand(strand_ratio)))

        print('\n--- MAPQ (mapped primary reads) ---')
        mean_mapq = float(np.mean(mapqs))
        high_qual = sum(1 for m in mapqs if m >= 30)
        print(f'Mean MAPQ: {mean_mapq:.1f}')
        print(f'MAPQ >= 30: {100*high_qual/len(mapqs):.1f}%')
        results.append(('Mean MAPQ', grade_min(mean_mapq, 40, 30)))

    print('\n--- Quality Summary ---')
    for label, grade in results:
        print(f'{label}: {grade}')
    fails = [label for label, grade in results if grade == 'FAIL']
    warns = [label for label, grade in results if grade == 'WARN']
    if fails:
        print('FAIL:', ', '.join(fails))
    if warns:
        print('WARN:', ', '.join(warns))
    if not fails and not warns:
        print('All metrics within normal range')
    return 1 if fails else 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Validate BAM alignment quality')
    parser.add_argument('bam', help='Input BAM file (index not required)')
    parser.add_argument('-n', '--sample-size', type=int, default=0,
                        help='Sample only the first N primary records (default 0 = whole file; see header note on bias)')
    args = parser.parse_args()

    sys.exit(validate_bam(args.bam, args.sample_size))
