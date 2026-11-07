#!/usr/bin/env python3
'''
Sashimi plot generation using ggsashimi.
Visualizes splicing events with read coverage and junction counts.

ggsashimi exits 0 when it silently drops missing BAMs, draws an empty region or hits an R error
(no figure written), so every call here is validated before and after the run.
'''
# Reference: ggsashimi 1.1.5 (GitHub guigolab/ggsashimi; needs R with ggplot2 < 3.5) | pysam 0.22+ | pandas 2.2+ | Verify API if version differs

import re
import subprocess
import numpy as np
import pandas as pd
import pysam
from pathlib import Path


def create_grouping_file(bam_files, conditions, output_file):
    '''
    Create a ggsashimi input TSV.

    ggsashimi expects column 1 = sample id/label, column 2 = BAM path,
    column 3 = group (used by -O for the overlay track and -C for the colour).

    Args:
        bam_files: List of BAM file paths
        conditions: List of group labels (same length as bam_files)
        output_file: Output TSV path
    '''
    groups = pd.DataFrame({
        'sample_id': [Path(b).stem for b in bam_files],
        'bam': bam_files,
        'group': conditions
    })
    groups.to_csv(output_file, sep='\t', index=False, header=False)
    return output_file


def write_palette(colors, output_file):
    '''One colour per line, assigned to groups in order of first appearance in the grouping TSV.'''
    Path(output_file).write_text('\n'.join(colors) + '\n')
    return output_file


def match_contig(bam_files, chrom):
    '''
    Return the contig name as spelled in the BAM headers (rMATS writes chrX, Ensembl BAMs use X).
    Raises if the name (with or without a "chr" prefix) is not in every BAM header.
    '''
    bare = chrom[3:] if chrom.startswith('chr') else chrom
    for cand in (chrom, bare, 'chr' + bare):
        found = True
        for bam_file in bam_files:
            with pysam.AlignmentFile(bam_file) as bam:
                if cand not in bam.references:
                    found = False
                    break
        if found:
            return cand
    raise ValueError(f'contig {chrom!r} is not in the header of every BAM')


def junction_counts(bam_file, chrom, start, end):
    '''
    Independent per-read intron counts {(donor, acceptor): n} inside chrom:start-end (1-based,
    donor = first intron base) - the same keys and counts ggsashimi draws before -M and -A.
    '''
    counts = {}
    with pysam.AlignmentFile(bam_file) as bam:
        for read in bam.fetch(chrom, start - 1, end):
            pos = read.reference_start + 1
            for op, length in read.cigartuples or []:
                if op == 3 and pos > start - 1 and pos + length < end:
                    counts[(pos, pos + length)] = counts.get((pos, pos + length), 0) + 1
                if op in (0, 2, 3, 7, 8):
                    pos += length
    return counts


def plot_sashimi(grouping_file, region, output_prefix, gtf_file, options=None):
    '''
    Generate sashimi plot for a genomic region and verify that the figure was written.

    Args:
        grouping_file: ggsashimi TSV (sample_id, bam, group columns)
        region: Genomic coordinates (chr:start-end, 1-based)
        output_prefix: Output file prefix
        gtf_file: GTF annotation file
        options: Dict of additional options: min_junc (-M, default 1; applied per sample
                 BEFORE aggregation), alpha, height, width, format, shrink, fix_y_scale,
                 overlay (default True: group column -> -O 3 -C 3), aggregate (mean_j default,
                 needs overlay), palette (file from write_palette), strand (-s; then two
                 files, <prefix>_+ and <prefix>_-)

    Returns:
        List of figure paths (all exist and are non-empty)
    '''
    if options is None:
        options = {}
    fmt = options.get('format', 'pdf')
    overlay = options.get('overlay', True)
    strand = options.get('strand', 'NONE')

    rows = [line.rstrip('\n').split('\t') for line in open(grouping_file) if line.strip()]
    # ggsashimi resolves relative BAM paths against the grouping file's directory
    bams = [r[1] if Path(r[1]).is_absolute() else str(Path(grouping_file).resolve().parent / r[1]) for r in rows]
    missing = [b for b in bams if not Path(b).is_file()]
    if missing:
        raise FileNotFoundError(f'BAM(s) missing (ggsashimi would drop them silently): {missing}')

    m = re.fullmatch(r'(.+):(-?[\d,]+)-([\d,]+)', region)  # start may be negative when an event sits near a contig start
    if m is None:
        raise ValueError(f'region must look like chr:start-end, got {region!r}')
    chrom, start, end = m.group(1), int(m.group(2).replace(',', '')), int(m.group(3).replace(',', ''))
    chrom = match_contig(bams, chrom)
    start = max(1, start)
    region = f'{chrom}:{start}-{end}'

    n_reads = 0
    for b in bams:
        with pysam.AlignmentFile(b) as bam:
            n_reads += bam.count(chrom, start - 1, end)
    if n_reads == 0:
        raise ValueError(f'no reads in {region} in any BAM (ggsashimi would write an empty figure)')
    min_junc = options.get('min_junc', 1)
    best = max((n for b in bams for n in junction_counts(b, chrom, start, end).values()), default=0)
    shrink = options.get('shrink', True)
    if best < min_junc:
        # ggsashimi exits 0 and writes a coverage-only figure with no arcs; with --shrink it crashes instead
        print(f'WARNING: no junction has >= {min_junc} reads (best {best}) in {region}: the figure will have no arcs; lower min_junc')
        if shrink:
            print('WARNING: dropping --shrink (ggsashimi crashes on it when no junction passes -M)')
            shrink = False

    cmd = [
        'ggsashimi.py',
        '-b', grouping_file,
        '-c', region,
        '-o', output_prefix,
        '-M', str(min_junc),
        '--alpha', str(options.get('alpha', 0.25)),
        '--height', str(options.get('height', 3)),
        '--width', str(options.get('width', 8)),
        '-F', fmt
    ]
    if gtf_file:
        cmd.extend(['-g', gtf_file])
    if strand != 'NONE':
        cmd.extend(['-s', strand])
    if overlay:
        cmd.extend(['-O', '3', '-C', '3'])
        cmd.extend(['-A', options.get('aggregate', 'mean_j')])  # mean, median, mean_j, median_j
    if options.get('palette'):
        cmd.extend(['-P', options['palette']])
    if shrink:
        cmd.append('--shrink')
    if options.get('fix_y_scale', True):
        cmd.append('--fix-y-scale')

    subprocess.run(cmd, check=True)

    expected = [f'{output_prefix}.{fmt}'] if strand == 'NONE' else [f'{output_prefix}_+.{fmt}', f'{output_prefix}_-.{fmt}']
    empty = [f for f in expected if not Path(f).is_file() or Path(f).stat().st_size == 0]
    if empty:
        raise RuntimeError(f'ggsashimi exited 0 but wrote no figure: {empty} (R error above? ggplot2 must be < 3.5)')
    print(f'Sashimi plot saved: {", ".join(expected)}')
    return expected


def batch_plot_rmats_events(rmats_file, grouping_file, gtf_file, output_dir,
                            n_top=20, fdr_cutoff=0.05, dpsi_cutoff=0.1, flank=500, min_junc=1, palette=None):
    '''
    Generate sashimi plots for top differential splicing events from rMATS.

    Args:
        rmats_file: rMATS SE.MATS.JC.txt output
        grouping_file: ggsashimi grouping file
        gtf_file: GTF annotation
        output_dir: Directory for output plots
        n_top: Number of top events to plot
        fdr_cutoff: FDR significance threshold
        dpsi_cutoff: Minimum |deltaPSI|
        flank: Bases of context added each side of upstreamES..downstreamEE
        min_junc: ggsashimi -M (per sample, before aggregation); 1 = draw every junction
        palette: colour file from write_palette (without it -C 3 gives R's default red/green)

    Raises RuntimeError after the loop if any event failed to plot.
    '''
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(rmats_file, sep='\t')

    # Filter significant events
    significant = df[
        (df['FDR'] < fdr_cutoff) &
        (df['IncLevelDifference'].abs() > dpsi_cutoff)
    ].copy()

    # Sort by significance and effect size
    significant['score'] = -significant['FDR'].apply(lambda x: max(x, 1e-300)).apply(np.log10) * significant['IncLevelDifference'].abs()
    significant = significant.nlargest(n_top, 'score')

    print(f'Plotting top {len(significant)} significant events')

    failed = []
    for _, event in significant.iterrows():
        chrom = event['chr']
        gene = event['geneSymbol']

        # Flanking exons plus context; plot_sashimi clamps the start to 1 and maps chrX <-> X
        region = f'{chrom}:{event["upstreamES"] - flank}-{event["downstreamEE"] + flank}'

        safe_gene = re.sub(r'[^A-Za-z0-9._-]', '_', str(gene))
        output_prefix = f'{output_dir}/{safe_gene}_{chrom}_{event["exonStart_0base"]}_{event["ID"]}'

        try:
            plot_sashimi(
                grouping_file, region, output_prefix, gtf_file,
                options={'shrink': True, 'fix_y_scale': True, 'min_junc': min_junc, 'palette': palette}
            )
        except (subprocess.CalledProcessError, ValueError, RuntimeError, FileNotFoundError) as e:
            print(f'Failed to plot {gene}: {e}')
            failed.append(gene)

    if failed:
        raise RuntimeError(f'{len(failed)} of {len(significant)} events failed: {failed}')


def plot_specific_event(grouping_file, gtf_file, chrom, start, end,
                        output_prefix, gene_name=None, flank=500, min_junc=1):
    '''
    Plot a specific genomic region with optional flanking sequence.

    Args:
        grouping_file: ggsashimi grouping file
        gtf_file: GTF annotation
        chrom: Chromosome
        start: Start position
        end: End position
        output_prefix: Output file prefix
        gene_name: Optional gene name for labeling
        flank: Base pairs to add on each side
        min_junc: ggsashimi -M (per sample, before aggregation); 1 = draw every junction
    '''
    region = f'{chrom}:{start - flank}-{end + flank}'

    return plot_sashimi(
        grouping_file, region, output_prefix, gtf_file,
        options={
            'shrink': True,
            'fix_y_scale': True,
            'min_junc': min_junc,
            'aggregate': 'mean_j',
            'height': 4,
            'width': 10
        }
    )


if __name__ == '__main__':
    # Example: Create grouping file and palette (one colour per group)
    bams = ['ctrl1.bam', 'ctrl2.bam', 'treat1.bam', 'treat2.bam']
    conditions = ['Control', 'Control', 'Treatment', 'Treatment']

    # create_grouping_file(bams, conditions, 'sashimi_groups.tsv')
    # write_palette(['#1f77b4', '#ff7f0e'], 'palette.txt')

    # Example: Plot single region
    # plot_sashimi('sashimi_groups.tsv', 'chr1:1000000-1010000', 'example', 'annotation.gtf',
    #              options={'palette': 'palette.txt'})

    # Example: Batch plot rMATS results
    # batch_plot_rmats_events(
    #     'rmats_output/SE.MATS.JC.txt',
    #     'sashimi_groups.tsv',
    #     'annotation.gtf',
    #     'sashimi_plots/'
    # )

    print('Configure BAM files and run functions to generate sashimi plots')
