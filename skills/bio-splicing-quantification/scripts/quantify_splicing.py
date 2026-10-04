#!/usr/bin/env python3
'''
Splicing quantification using SUPPA2 and rMATS-turbo.
Demonstrates PSI calculation from RNA-seq data.
'''
# Reference: kallisto 0.50+, pandas 2.2+ | Verify API if version differs

import subprocess
import pandas as pd
from pathlib import Path


def run_suppa2_quantification(gtf_file, tpm_file, output_prefix, event_types=None):
    '''
    Quantify splicing using SUPPA2 from transcript TPM.

    Args:
        gtf_file: Path to GTF annotation
        tpm_file: Tab-separated TPM file (transcripts x samples); header has sample
            names only, see write_suppa_tpm
        output_prefix: Prefix for output files
        event_types: List of event types (default: all)
    '''
    if event_types is None:
        event_types = ['SE', 'SS', 'MX', 'RI', 'FL']

    # Step 1: Generate splicing events from annotation
    subprocess.run([
        'suppa.py', 'generateEvents',
        '-i', gtf_file,
        '-o', output_prefix,
        '-f', 'ioe',
        '-e'] + event_types,
        check=True
    )

    # Step 2: Calculate PSI for each output event code.
    # generateEvents writes one file per OUTPUT code: SS -> A5 and A3; FL -> AF and AL.
    event_codes = {'SE': ['SE'], 'SS': ['A5', 'A3'], 'MX': ['MX'], 'RI': ['RI'], 'FL': ['AF', 'AL']}
    psi_files = {}

    for et in event_types:
        for code in event_codes.get(et, [et]):
            ioe_file = f'{output_prefix}_{code}_strict.ioe'
            psi_output = f'{output_prefix}_psi_{code}'

            if Path(ioe_file).exists():
                subprocess.run([
                    'suppa.py', 'psiPerEvent',
                    '-i', ioe_file,
                    '-e', tpm_file,
                    '-o', psi_output
                ], check=True)
                psi_files[code] = f'{psi_output}.psi'

    return psi_files


def write_suppa_tpm(tpm, path):
    '''
    Write a transcript x sample TPM DataFrame in the layout SUPPA2 requires.

    The header line holds only the sample names (no cell above the transcript
    IDs). A pandas default header (index name in the first cell) makes
    psiPerEvent fail with "6 expected, 5 given" and write nothing.
    '''
    with open(path, 'w', newline='\n') as handle:
        handle.write('\t'.join(str(c) for c in tpm.columns) + '\n')
        tpm.to_csv(handle, sep='\t', header=False, lineterminator='\n')


def filter_reliable_events(psi_file, min_samples_with_coverage=0.5, psi_range=(0.05, 0.95)):
    '''
    Filter PSI matrix for reliable, non-constitutive events.

    Args:
        psi_file: Path to PSI file from SUPPA2
        min_samples_with_coverage: Minimum fraction of samples with non-NA PSI
        psi_range: (low, high); keep events whose mean PSI lies strictly inside.
            This is a filter applied here, not a SUPPA2 or rMATS default.
    '''
    psi = pd.read_csv(psi_file, sep='\t', index_col=0)

    # NA indicates insufficient reads (TPM) for PSI calculation
    coverage_frac = psi.notna().mean(axis=1)
    reliable = psi[coverage_frac >= min_samples_with_coverage]

    # Remove near-constitutive events (mean PSI close to 0 or 1)
    low, high = psi_range
    mean_psi = reliable.mean(axis=1)
    variable = reliable[(mean_psi > low) & (mean_psi < high)]

    print(f'Total events: {len(psi)}')
    print(f'After coverage filter: {len(reliable)}')
    print(f'After variability filter: {len(variable)}')

    return variable


# Event-specific coordinate columns in rMATS *.MATS.JC.txt / *.MATS.JCEC.txt.
RMATS_COORD_COLUMNS = {
    'SE': ['exonStart_0base', 'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'],
    'A5SS': ['longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE'],
    'A3SS': ['longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE'],
    'MXE': ['1stExonStart_0base', '1stExonEnd', '2ndExonStart_0base', '2ndExonEnd',
            'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'],
    'RI': ['riExonStart_0base', 'riExonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'],
}


def _csv_values(cell, cast=float):
    '''Split an rMATS comma-separated per-replicate cell; NA entries are dropped.'''
    return [cast(x) for x in str(cell).split(',') if x not in ('NA', 'nan', '')]


def parse_rmats_output(rmats_dir, event_type='SE', min_junction_reads=20, counts='JC'):
    '''
    Parse rMATS output and calculate mean PSI for any of the five event types.

    Args:
        rmats_dir: Directory containing rMATS output
        event_type: SE, A5SS, A3SS, MXE, or RI
        min_junction_reads: Minimum IJC+SJC in EVERY replicate of BOTH groups
        counts: 'JC' (junction reads only) or 'JCEC' (adds exon-body reads)

    mean_PSI pools IncLevel1 and IncLevel2 only. IncLevelDifference is a
    group contrast, not a PSI, and is never included.
    '''
    if event_type not in RMATS_COORD_COLUMNS:
        raise ValueError(f'event_type must be one of {sorted(RMATS_COORD_COLUMNS)}')
    df = pd.read_csv(Path(rmats_dir) / f'{event_type}.MATS.{counts}.txt', sep='\t')

    def row_stats(row):
        inc1, inc2 = _csv_values(row['IncLevel1']), _csv_values(row['IncLevel2'])
        pooled = inc1 + inc2
        reads = [int(i) + int(s)
                 for ij, sj in (('IJC_SAMPLE_1', 'SJC_SAMPLE_1'), ('IJC_SAMPLE_2', 'SJC_SAMPLE_2'))
                 for i, s in zip(str(row[ij]).split(','), str(row[sj]).split(','))]
        return pd.Series({
            'mean_PSI_group1': sum(inc1) / len(inc1) if inc1 else None,
            'mean_PSI_group2': sum(inc2) / len(inc2) if inc2 else None,
            'mean_PSI': sum(pooled) / len(pooled) if pooled else None,
            'min_reads_per_replicate': min(reads),
        })

    df = pd.concat([df, df.apply(row_stats, axis=1)], axis=1)
    reliable = df[df['min_reads_per_replicate'] >= min_junction_reads].copy()

    print(f'Total {event_type} events: {len(df)}')
    print(f'Reliable events (>={min_junction_reads} reads in every replicate): {len(reliable)}')

    return reliable[['ID', 'GeneID', 'geneSymbol', 'chr', 'strand'] + RMATS_COORD_COLUMNS[event_type]
                    + ['mean_PSI_group1', 'mean_PSI_group2', 'mean_PSI', 'min_reads_per_replicate']]


if __name__ == '__main__':
    # Example usage with SUPPA2
    # Requires: annotation.gtf, transcript_tpm.tsv from Salmon/kallisto
    gtf = 'annotation.gtf'
    tpm = 'transcript_tpm.tsv'

    # Run SUPPA2 quantification
    # psi_files = run_suppa2_quantification(gtf, tpm, 'splicing_events')

    # Filter for reliable events
    # reliable_se = filter_reliable_events('splicing_events_psi_SE.psi')  # psi_range=(0.05, 0.95)
    # print(reliable_se.head())

    # Example with rMATS output
    # reliable_se = parse_rmats_output('rmats_output/', 'SE', min_junction_reads=20)
    # print(reliable_se.head())

    print('Run with actual data files to quantify splicing events')
