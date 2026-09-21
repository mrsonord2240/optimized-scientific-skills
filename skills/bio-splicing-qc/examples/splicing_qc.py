#!/usr/bin/env python3
'''
Quality control for splicing analysis: RSeQC junction annotation / saturation,
fragment-level junction counts with the real anchor (overhang) of each junction,
and MaxEntScan splice-site scoring.

Checked 2026-09-20 on RSeQC 5.0.5, pysam 0.24.1, pandas 2.3.3, maxentpy 0.0.2
(see test_splicing_qc.py).

CLI:
  python splicing_qc.py annotation sample.bam genes.bed12 out_prefix
  python splicing_qc.py saturation sample.bam genes.bed12 out_prefix
  python splicing_qc.py junctions  sample.bam [--min-overhang 8]
  python splicing_qc.py report     sample.bam genes.bed12 out_prefix
  python splicing_qc.py sites --donor CAGGTAAGT --acceptor TTTTTTTTTTTTTTCCTTAGGAG

RSeQC needs a BED12 gene model whose contig names match the BAM. The plotting step needs
Rscript; every RSeQC call here passes --skip-plot unless plot=True and Rscript is on PATH.
'''

import argparse
import math
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam

ALIGNED_OPS = (0, 7, 8)          # M, =, X
REF_OPS = (0, 2, 3, 7, 8)        # M, D, N, =, X consume the reference


class SplicingQCError(RuntimeError):
    '''Input problem that RSeQC would otherwise hide behind exit code 0 or 1.'''


class NoSplicedReadsError(SplicingQCError):
    '''The BAM has no spliced (N) reads passing the filters.'''


def check_contigs(bam_path, bed_path, need_bed12=True):
    '''
    Raise if BAM and BED share no contig name, or if the gene model is not BED12.
    RSeQC exits 0 on both: a chr1 vs 1 mismatch gives "Total 0 usable reads" and flat zero
    saturation curves; a BED6 gene model classifies every junction as complete_novel.
    '''
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as bam:
        bam_contigs = set(bam.references)
    bed_contigs, narrow = set(), 0
    with open(bed_path) as fh:
        for line in fh:
            if line.strip() and not line.startswith(('track', 'browser', '#')):
                fields = line.rstrip('\r\n').split('\t')
                bed_contigs.add(fields[0])
                narrow += len(fields) < 12
    if need_bed12 and narrow:
        raise SplicingQCError(f'{bed_path}: {narrow} lines have fewer than 12 columns; junction_annotation and '
                              'junction_saturation need a BED12 gene model (see SKILL.md, "Reference files")')
    shared = bam_contigs & bed_contigs
    if not shared:
        raise SplicingQCError(
            f'no contig name shared by BAM ({sorted(bam_contigs)[:3]}...) and BED '
            f'({sorted(bed_contigs)[:3]}...): rename one (chr1 vs 1) before running RSeQC')
    if len(bed_contigs - bam_contigs) > len(bed_contigs) / 2:
        print(f'warning: {len(bed_contigs - bam_contigs)} of {len(bed_contigs)} BED contigs are absent from the BAM',
              file=sys.stderr)
    return shared


def _run_rseqc(tool, bam_file, bed_file, prefix, plot, extra=()):
    cmd = [tool, '-i', str(bam_file), '-r', str(bed_file), '-o', str(prefix), *extra]
    if not (plot and shutil.which('Rscript')):
        cmd.append('--skip-plot')   # without Rscript, plotting makes RSeQC exit 1 after writing its tables
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SplicingQCError(f'{tool} failed (exit {proc.returncode}): {proc.stderr.strip()[-400:]}')
    return proc


def junction_annotation(bam_file, bed_file, output_prefix, plot=False):
    '''
    Run RSeQC junction_annotation.py and return known / partial-novel / complete-novel fractions.

    RSeQC calls a junction "annotated" when its donor AND acceptor are each in the gene model, so a
    skipping junction between two annotated exons counts as known. Both denominators are returned:
    'read_*' (reads, what the SKILL.md thresholds use) and 'junction_*' (distinct junctions).
    '''
    check_contigs(bam_file, bed_file)
    _run_rseqc('junction_annotation.py', bam_file, bed_file, output_prefix, plot)
    xls = Path(f'{output_prefix}.junction.xls')
    try:
        junc = pd.read_csv(xls, sep='\t')      # RSeQC writes a 0-byte file when there are no spliced reads
    except (FileNotFoundError, pd.errors.EmptyDataError):
        junc = pd.DataFrame()
    if junc.empty:
        raise NoSplicedReadsError(f'{bam_file}: RSeQC found no spliced reads at MAPQ >= 30 (unspliced or DNA-like library, '
                                  'or -q too strict for this aligner)')
    junc['annotation'] = junc['annotation'].str.strip()      # RSeQC writes ' annotated' with a leading space
    reads = junc.groupby('annotation')['read_count'].sum()
    njunc = junc.groupby('annotation').size()
    novel_classes = ('partial_novel', 'complete_novel')
    out = {'total_reads': int(reads.sum()), 'total_junctions': int(njunc.sum())}
    for label, series in (('read', reads), ('junction', njunc)):
        total = series.sum()
        known = series.get('annotated', 0) / total
        novel = sum(series.get(c, 0) for c in novel_classes) / total
        assert abs(known + novel - 1) < 1e-9, f'unexpected annotation classes: {list(series.index)}'
        out[f'{label}_known'] = float(known)
        out[f'{label}_novel'] = float(novel)
    out['status'] = ('Healthy' if out['read_known'] >= 0.8 else
                     'Acceptable' if out['read_known'] >= 0.6 else 'Suspect or interesting')
    return out


def read_saturation_curve(r_file):
    '''Parse x (percent of reads), y (known), z (all), w (novel) from RSeQC's *.junctionSaturation_plot.r.'''
    text = Path(r_file).read_text()
    curve = {}
    for name in 'xyzw':
        m = re.search(rf'^{name}=c\(([^)]*)\)', text, re.M)
        if not m:
            raise SplicingQCError(f'{r_file}: no {name}=c(...) vector; not a RSeQC saturation script')
        curve[name] = [float(v) for v in m.group(1).split(',')]
    return {'percent': curve['x'], 'known': curve['y'], 'all': curve['z'], 'novel': curve['w']}


def junction_saturation(bam_file, bed_file, output_prefix, plot=False, step=5, lo=5, hi=100):
    '''
    Run RSeQC junction_saturation.py and apply the plateau rule (growth from 80% to 100% of reads < 2%).

    RSeQC shuffles splice events without a seed: repeat runs differ by a few junctions.
    Returns the curves, growth of the known and all-junction curves, and the verdict.
    '''
    if lo != step:
        raise ValueError('RSeQC samples correctly only when lo == step (-l equals -s); a finer curve means a smaller step, e.g. lo=step=2')
    check_contigs(bam_file, bed_file)
    _run_rseqc('junction_saturation.py', bam_file, bed_file, output_prefix, plot,
               extra=('-l', str(lo), '-u', str(hi), '-s', str(step)))
    curve = read_saturation_curve(f'{output_prefix}.junctionSaturation_plot.r')
    if not any(curve['all']) or curve['all'][-1] == 0:
        raise NoSplicedReadsError(f'{bam_file}: saturation curve is all zeros (no spliced reads, or BED/BAM mismatch)')
    pct = curve['percent']
    i80 = max(i for i, p in enumerate(pct) if p <= 80) if pct[0] <= 80 else 0
    growth = {}
    for name in ('known', 'all'):
        base = curve[name][i80]
        growth[name] = (curve[name][-1] - base) / base if base else float('nan')
    curve['growth_80_100'] = growth
    curve['verdict'] = ('PLATEAU' if growth['known'] < 0.02 else 'STILL RISING (under-sequenced)')
    return curve


def junction_stats(bam_path, min_overhang=8, min_mapq=30, unique_only=True):
    '''
    Fragment-level junction support from a BAM (works without an index).

    For every N operation the overhang is the aligned length (M/=/X) of the block immediately
    left and right of it, so a read 30M1000N4M800N66M has overhangs 4 (both junctions), not 30.
    Secondary and supplementary records are skipped; unique_only keeps NH==1 (or MAPQ >= min_mapq when
    NH is absent); the two mates of a pair count once per junction. Reads with overhang < min_overhang
    on either side are counted in 'reads_all' but not in 'reads'.

    Returns {(contig, intron_start_0based, intron_end): {'reads', 'reads_all', 'min_overhang'}}.
    Memory grows with the number of spliced read names; for very deep BAMs use STAR SJ.out.tab or regtools.
    '''
    stats = defaultdict(lambda: {'reads': 0, 'reads_all': 0, 'min_overhang': None})
    seen_all, seen_ok = set(), set()
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.cigartuples:
                continue
            if unique_only:
                if read.has_tag('NH'):
                    if read.get_tag('NH') != 1:
                        continue
                elif read.mapping_quality < min_mapq:
                    continue
            # split the CIGAR into aligned blocks separated by N
            blocks, introns, cur, pos = [], [], 0, read.reference_start
            for op, length in read.cigartuples:
                if op == 3:
                    blocks.append(cur)
                    introns.append((pos, pos + length))
                    cur = 0
                elif op in ALIGNED_OPS:
                    cur += length
                if op in REF_OPS:
                    pos += length
            blocks.append(cur)
            for k, (start, end) in enumerate(introns):
                key = (read.reference_name, start, end)
                overhang = min(blocks[k], blocks[k + 1])
                frag = (read.query_name, key) if read.is_paired else None
                rec = stats[key]
                if frag is None or frag not in seen_all:
                    rec['reads_all'] += 1
                    if frag is not None:
                        seen_all.add(frag)
                if overhang >= min_overhang and (frag is None or frag not in seen_ok):
                    rec['reads'] += 1
                    if frag is not None:
                        seen_ok.add(frag)
                rec['min_overhang'] = overhang if rec['min_overhang'] is None else min(rec['min_overhang'], overhang)
    return dict(stats)


def summarize_junctions(stats, min_reads=10):
    '''Counts of junctions supported by >= min_reads reads that pass the overhang filter.'''
    total = len(stats)
    supported = [k for k, v in stats.items() if v['reads'] > 0]
    ge = sum(1 for v in stats.values() if v['reads'] >= min_reads)
    return {
        'total_junctions': total,
        'junctions_with_anchored_reads': len(supported),
        f'junctions_ge_{min_reads}_reads': ge,
        'pct_ge_min_reads': 100 * ge / total if total else 0.0,
        'reads_total': sum(v['reads_all'] for v in stats.values()),
        'reads_anchored': sum(v['reads'] for v in stats.values()),
    }


def score_splice_sites(sequences_5ss, sequences_3ss):
    '''
    MaxEntScan scores, one output per input (NaN for invalid input, so indices stay aligned).

    5'ss: 9 nt (3 exon + 6 intron); 3'ss: 23 nt (20 intron + 3 exon, ending ...AG|exon).
    maxentpy exits the process on a wrong length (SystemExit) and raises KeyError on N/U, so both are
    checked here. Requires: pip install maxentpy
    '''
    from maxentpy.maxent import score5, score3

    def one(fn, seq, length):
        seq = seq.upper()
        if len(seq) != length or set(seq) - set('ACGT'):
            return math.nan
        try:
            return float(fn(seq))
        except (SystemExit, KeyError, ValueError):
            return math.nan

    return ([one(score5, s, 9) for s in sequences_5ss],
            [one(score3, s, 23) for s in sequences_3ss])


def generate_qc_report(bam_file, bed_file, output_prefix, plot=False):
    '''Annotation, saturation and junction-support QC for one BAM.'''
    print(f'Analyzing: {bam_file}')
    print('=' * 50)

    print('\n1. Known vs novel junctions (RSeQC junction_annotation)')
    ann = junction_annotation(bam_file, bed_file, output_prefix, plot)
    print(f"  known {ann['read_known']:.1%} of {ann['total_reads']} reads "
          f"({ann['junction_known']:.1%} of {ann['total_junctions']} junctions) -> {ann['status']}")

    print('\n2. Junction saturation (RSeQC junction_saturation)')
    sat = junction_saturation(bam_file, bed_file, output_prefix, plot)
    print(f"  known curve growth 80->100% of reads: {sat['growth_80_100']['known']:.1%} -> {sat['verdict']}")

    print('\n3. Junction read support (pysam, fragment-level, unique reads)')
    summary = summarize_junctions(junction_stats(bam_file))
    for key, value in summary.items():
        print(f'  {key}: {value:.1f}' if isinstance(value, float) else f'  {key}: {value}')
    pct = summary['pct_ge_min_reads']
    print('  Junction coverage:', 'GOOD (>=50% of junctions have >=10 reads)' if pct >= 50 else
          'ACCEPTABLE (30-50%)' if pct >= 30 else 'POOR (<30%): consider deeper sequencing')
    return {'annotation': ann, 'saturation': sat, 'support': summary}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='cmd', required=True)
    for name in ('annotation', 'saturation', 'report'):
        p = sub.add_parser(name)
        p.add_argument('bam')
        p.add_argument('bed12')
        p.add_argument('prefix')
        p.add_argument('--plot', action='store_true', help='draw the RSeQC PDFs (needs Rscript)')
    p = sub.add_parser('junctions')
    p.add_argument('bam')
    p.add_argument('--min-overhang', type=int, default=8)
    p = sub.add_parser('sites')
    p.add_argument('--donor', action='append', default=[])
    p.add_argument('--acceptor', action='append', default=[])
    args = parser.parse_args(argv)

    if args.cmd == 'annotation':
        print(junction_annotation(args.bam, args.bed12, args.prefix, args.plot))
    elif args.cmd == 'saturation':
        print(junction_saturation(args.bam, args.bed12, args.prefix, args.plot))
    elif args.cmd == 'report':
        generate_qc_report(args.bam, args.bed12, args.prefix, args.plot)
    elif args.cmd == 'junctions':
        print(summarize_junctions(junction_stats(args.bam, min_overhang=args.min_overhang)))
    else:
        s5, s3 = score_splice_sites(args.donor, args.acceptor)
        print("5'ss:", [round(v, 2) for v in s5], "3'ss:", [round(v, 2) for v in s3])


if __name__ == '__main__':
    try:
        main()
    except SplicingQCError as exc:
        sys.exit(f'error: {exc}')
