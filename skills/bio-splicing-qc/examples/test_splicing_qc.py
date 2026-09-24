#!/usr/bin/env python3
'''
Self-contained test for splicing_qc.py: builds a tiny BAM and BED12 with planted truth, then asserts.
Run:  python test_splicing_qc.py
The RSeQC checks are skipped when junction_annotation.py is not on PATH, the MaxEntScan checks
when maxentpy is not importable. Needs pysam and pandas.
'''
import re
import shutil
import sys
import tempfile
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).parent))
import splicing_qc as sq  # noqa: E402

CHROM, CHROM_LEN = 'chrT', 20000
# gene model: exons [1000,1200) [2000,2200) [3000,3200) [4000,4200) on +
EXONS = [(1000, 1200), (2000, 2200), (3000, 3200), (4000, 4200)]


def write_bed12(path, chrom=CHROM):
    starts = [e[0] - EXONS[0][0] for e in EXONS]
    sizes = [e[1] - e[0] for e in EXONS]
    path.write_text('\t'.join([chrom, str(EXONS[0][0]), str(EXONS[-1][1]), 'geneT', '0', '+',
                               str(EXONS[0][0]), str(EXONS[-1][1]), '0', str(len(EXONS)),
                               ','.join(map(str, sizes)) + ',', ','.join(map(str, starts)) + ',']) + '\n')


def make_read(header, name, pos, cigar, mapq=60, nh=1, flag=0, mate=None):
    r = pysam.AlignedSegment(header)
    r.query_name, r.reference_id, r.reference_start = name, 0, pos
    r.mapping_quality, r.flag = mapq, flag
    r.cigarstring = cigar
    qlen = sum(int(n) for n, op in re.findall(r'(\d+)([MIS=X])', cigar))
    r.query_sequence = 'A' * qlen
    r.query_qualities = pysam.qualitystring_to_array('I' * qlen)
    r.set_tag('NH', nh)
    if mate:
        r.next_reference_id, r.next_reference_start = 0, mate
    return r


def write_bam(path, reads, header, index=True):
    reads = sorted(reads, key=lambda r: r.reference_start)
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for r in reads:
            out.write(r)
    if index:
        pysam.index(str(path))


def main():
    tmp = Path(tempfile.mkdtemp(prefix='splicing_qc_test_'))
    header = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6', 'SO': 'coordinate'},
                                              'SQ': [{'SN': CHROM, 'LN': CHROM_LEN}]})
    bed = tmp / 'genes.bed12'
    write_bed12(bed)

    # planted spliced reads (single-end, 100 nt, junction in the middle unless stated)
    reads = []
    def add(n, pos, cigar, prefix, **kw):
        for i in range(n):
            reads.append(make_read(header, f'{prefix}{i}', pos, cigar, **kw))
    add(60, 1150, '50M800N50M', 'kn1')      # 1200-2000, known
    add(20, 2150, '50M800N50M', 'kn2')      # 2200-3000, known
    add(10, 1150, '50M1300N50M', 'part')    # 1200-2500: donor known, acceptor novel (partial)
    add(10, 1450, '50M1050N50M', 'nov')     # 1500-2550: both novel
    add(6, 1150, '50M800N50M', 'mm', nh=4, mapq=1)   # multi-mapper on a known junction: must be ignored
    spliced_bam = tmp / 'spliced.bam'
    write_bam(spliced_bam, reads, header)

    # 1. overhang from the adjacent block only, minimum-overhang filter, =/X ops
    micro = [make_read(header, f'micro{i}', 1000, '30M1000N4M800N66M') for i in range(12)]
    micro += [make_read(header, f'eqx{i}', 5000, '20=1X29=500N50=') for i in range(3)]
    micro_bam = tmp / 'micro.bam'
    write_bam(micro_bam, micro, header, index=False)          # unindexed on purpose
    st = sq.junction_stats(str(micro_bam), min_overhang=8)
    j1, j2 = (CHROM, 1030, 2030), (CHROM, 2034, 2834)
    assert st[j1]['min_overhang'] == 4 and st[j2]['min_overhang'] == 4, st
    assert st[j1]['reads_all'] == 12 and st[j1]['reads'] == 0, 'overhang 4 < 8 must not count'
    assert sq.junction_stats(str(micro_bam), min_overhang=4)[j1]['reads'] == 12
    assert (CHROM, 5050, 5550) in st and st[(CHROM, 5050, 5550)]['reads'] == 3, 'junction start must count =/X ops'
    print('ok  overhang, min_overhang and =/X coordinates (unindexed BAM)')

    # 2. multi-mappers ignored; the two mates of a pair count once
    pair = []
    for i in range(5):
        pair.append(make_read(header, f'p{i}', 1150, '50M800N50M', flag=1 | 2 | 64 | 32, mate=1180))
        pair.append(make_read(header, f'p{i}', 1180, '20M800N80M', flag=1 | 2 | 128 | 16, mate=1150))
    pair_bam = tmp / 'pairs.bam'
    write_bam(pair_bam, pair, header)
    ps = sq.junction_stats(str(pair_bam))
    assert ps[(CHROM, 1200, 2000)]['reads'] == 5, ps        # 10 records, 5 fragments
    ss = sq.junction_stats(str(spliced_bam))
    assert ss[(CHROM, 1200, 2000)]['reads'] == 60, 'NH=4 reads must be excluded'
    summary = sq.summarize_junctions(ss, min_reads=10)
    assert summary['total_junctions'] == 4 and summary['junctions_ge_10_reads'] == 4, summary
    print('ok  fragment-level counting, NH filter')

    # 3. BAM without spliced reads: clear result, no crash
    plain = [make_read(header, f'u{i}', 1000, '100M') for i in range(20)]
    plain_bam = tmp / 'plain.bam'
    write_bam(plain_bam, plain, header)
    assert sq.junction_stats(str(plain_bam)) == {}
    assert sq.summarize_junctions({})['pct_ge_min_reads'] == 0.0
    print('ok  BAM without spliced reads')

    # 4. class-level support joins RSeQC coordinates to the helper output without dropping zero support
    xls = tmp / 'classes.junction.xls'
    xls.write_text('chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation\n'
                   f'{CHROM}\t1200\t2000\t60\t annotated\n'
                   f'{CHROM}\t1200\t2500\t10\t partial_novel\n'
                   f'{CHROM}\t1500\t2550\t10\t complete_novel\n')
    cls = sq.junction_class_support(xls, ss).set_index('annotation')
    assert cls.loc['annotated', 'junctions'] == 1
    assert cls.loc['annotated', 'median_anchored_fragments'] == 60
    assert cls.loc['partial_novel', 'median_anchored_fragments'] == 10
    assert cls.loc['complete_novel', 'median_anchored_fragments'] == 10
    print('ok  annotation-class support diagnostic')

    # 5. splice-site scoring keeps one output per input
    try:
        import maxentpy  # noqa: F401
    except ImportError:
        print('skip MaxEntScan (maxentpy not installed)')
    else:
        s5, s3 = sq.score_splice_sites(['CAGGTAAGT', 'CAGGTAA', 'NAGGTAAGT', 'CAGATAAGT'],
                                       ['TTTTTTTTTTTTTTCCTTAGGAG', 'TTTTTTTTTTTTTTTTTTTTCAG'])
        assert len(s5) == 4 and len(s3) == 2
        assert s5[0] > 8 > 5 > s5[3] and s5[1] != s5[1] and s5[2] != s5[2], s5    # NaN for invalid input
        assert s3[0] > 8, s3                     # real acceptor scores ~11.6
        assert s3[1] < 0, s3                     # no AG in the intron end: negative
        print('ok  MaxEntScan scoring aligned with inputs', [round(v, 2) for v in s5], [round(v, 2) for v in s3])

    # 5. RSeQC (skipped when not installed)
    if not shutil.which('junction_annotation.py'):
        print('skip RSeQC checks (junction_annotation.py not on PATH)')
    else:
        ann = sq.junction_annotation(str(spliced_bam), str(bed), str(tmp / 'ann'))       # no Rscript needed
        # known: 80 reads, partial 10, complete novel 10 (multi-mappers filtered by MAPQ)
        assert ann['total_reads'] == 100 and abs(ann['read_known'] - 0.8) < 1e-9, ann
        assert abs(ann['read_known'] + ann['read_novel'] - 1) < 1e-9
        assert ann['status'] == 'Healthy' and abs(ann['junction_known'] - 0.5) < 1e-9, ann
        try:
            sq.junction_annotation(str(plain_bam), str(bed), str(tmp / 'ann_plain'))
        except sq.NoSplicedReadsError:
            pass
        else:
            raise AssertionError('expected NoSplicedReadsError for a BAM without spliced reads')
        bad_bed = tmp / 'genes_nochr.bed12'
        write_bed12(bad_bed, chrom='T')
        try:
            sq.junction_annotation(str(spliced_bam), str(bad_bed), str(tmp / 'ann_bad'))
        except sq.SplicingQCError:
            pass
        else:
            raise AssertionError('expected SplicingQCError for a contig-name mismatch')
        bed6 = tmp / 'genes.bed6'
        bed6.write_text('\t'.join([CHROM, '1000', '4200', 'geneT', '0', '+']) + '\n')
        try:
            sq.junction_annotation(str(spliced_bam), str(bed6), str(tmp / 'ann_bed6'))
        except sq.SplicingQCError:
            pass
        else:
            raise AssertionError('expected SplicingQCError for a BED6 gene model')
        print('ok  known/novel fractions, empty BAM, contig mismatch, BED6')

        sat = sq.junction_saturation(str(spliced_bam), str(bed), str(tmp / 'sat'))
        assert sat['all'][-1] == 4 and sat['known'][-1] == 2, sat
        assert sat['verdict'].startswith('PLATEAU') or sat['verdict'].startswith('STILL'), sat
        print('ok  saturation curve parsed', sat['growth_80_100'], sat['verdict'])
        try:
            sq.junction_saturation(str(spliced_bam), str(bed), str(tmp / 'sat_bad'), lo=80, step=5)
        except ValueError:
            pass
        else:
            raise AssertionError('expected ValueError for lo != step')
        print('ok  lo != step refused')

    shutil.rmtree(tmp)
    print('all checks passed')


if __name__ == '__main__':
    main()
