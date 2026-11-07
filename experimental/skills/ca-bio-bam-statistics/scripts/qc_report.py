#!/usr/bin/env python3
'''Generate QC report from BAM file'''
# Checked on: pysam 0.24.1, samtools 1.24 | Verify API if version differs
# Counts follow `samtools flagstat` (QC-passed column): rates use primary alignments only;
# secondary/supplementary records and QC-failed reads are counted separately.
# Usage: qc_report.py <input.bam|cram> [reference.fa]   (CRAM needs the reference FASTA)

import pysam
import sys

# Templates >= this are dropped from the insert summary (samtools stats -i 8000 counts them at 8000 instead,
# so the means differ on long-insert libraries); not a 1000 bp library limit
MAX_INSERT = 8000

def pct(part, whole):
    return f'{part / whole * 100:.2f}%' if whole else 'n/a'

def qc_report(bam_path, reference=None):
    s = {'records': 0, 'secondary': 0, 'supplementary': 0, 'primary': 0, 'qcfail': 0,
         'mapped': 0, 'paired': 0, 'proper_pair': 0, 'duplicate': 0}
    insert_sizes = []

    is_cram = bam_path.lower().endswith('.cram')
    try:
        # check_sq=False: unaligned BAM/uBAM has no @SQ lines; reference is only used for CRAM
        with pysam.AlignmentFile(bam_path, 'rb', check_sq=False, reference_filename=reference) as bam:
            for read in bam:
                s['records'] += 1
                if read.is_secondary or read.is_supplementary:
                    s['secondary'] += read.is_secondary
                    s['supplementary'] += read.is_supplementary
                    continue
                s['primary'] += 1
                if read.is_qcfail:
                    s['qcfail'] += 1
                    continue
                if not read.is_unmapped:
                    s['mapped'] += 1
                if read.is_paired:
                    s['paired'] += 1
                if read.is_proper_pair:
                    s['proper_pair'] += 1
                    if read.is_read1 and 0 < read.template_length < MAX_INSERT:
                        insert_sizes.append(read.template_length)
                if read.is_duplicate:
                    s['duplicate'] += 1
    except (OSError, ValueError, NotImplementedError) as e:
        hint = ''
        if is_cram:
            hint = (' (check that the reference FASTA is the one the CRAM was written against: see the htslib MD5 message above)'
                    if reference else ' (CRAM cannot be decoded without its reference: pass the FASTA as 2nd argument)')
        sys.exit(f'qc_report.py: cannot read {bam_path}: {e}{hint}')

    passed = s['primary'] - s['qcfail']
    print(f'=== QC Report: {bam_path} ===\n')
    print(f'Total records:     {s["records"]:,} (primary {s["primary"]:,} + secondary {s["secondary"]:,} '
          f'+ supplementary {s["supplementary"]:,})')
    print(f'QC-failed primary: {s["qcfail"]:,} (excluded below)')
    print(f'QC-passed primary: {passed:,}')
    if not passed:
        print('No QC-passed primary reads: nothing to report')
        return
    print(f'Mapped:            {s["mapped"]:,} ({pct(s["mapped"], passed)} of QC-passed primary)')
    print(f'Properly paired:   {s["proper_pair"]:,} ({pct(s["proper_pair"], s["paired"])} of primary paired reads)')
    print(f'Duplicates:        {s["duplicate"]:,} ({pct(s["duplicate"], passed)} of QC-passed primary)')

    if insert_sizes:
        mean_insert = sum(insert_sizes) / len(insert_sizes)
        insert_sizes.sort()
        median_insert = insert_sizes[len(insert_sizes) // 2]
        print(f'\nInsert size (mean): {mean_insert:.0f}')
        print(f'Insert size (median): {median_insert}')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: qc_report.py <input.bam|cram> [reference.fa]')
        sys.exit(1)

    qc_report(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
