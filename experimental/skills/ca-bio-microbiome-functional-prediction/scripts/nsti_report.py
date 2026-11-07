#!/usr/bin/env python3
"""Report the NSTI distribution and the fraction of reads dropped by the --max_nsti gate.

Reference: PICRUSt2 2.6.3, pandas 2.2+ | Verify API if version differs
Inputs:  PICRUSt2 output dir (contains combined_marker_predicted_and_nsti.tsv.gz),
         ASV table TSV (ASVs x samples, no '# Constructed from biom file' line),
         optional --max-nsti (default 2.0, the PICRUSt2 default).
Usage:   python scripts/nsti_report.py picrust2_out asv_table.tsv [--max-nsti 2.0]
"""
import argparse

import pandas as pd

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument('picrust2_out', help='PICRUSt2 output directory')
ap.add_argument('asv_table', help='ASV abundance TSV, ASVs x samples')
ap.add_argument('--max-nsti', type=float, default=2.0,
                help='PICRUSt2 default; ASVs above this are dropped before metagenome inference')
args = ap.parse_args()

nsti = pd.read_csv(f'{args.picrust2_out}/combined_marker_predicted_and_nsti.tsv.gz', sep='\t')   # cols: sequence, metadata_NSTI
asv_counts = pd.read_csv(args.asv_table, sep='\t', index_col=0)                # ASVs x samples
nsti = nsti.set_index('sequence')
reads_per_asv = asv_counts.sum(axis=1)

max_nsti = args.max_nsti
dropped = nsti.index[nsti['metadata_NSTI'] > max_nsti]
reads_dropped_frac = reads_per_asv.reindex(dropped).sum() / reads_per_asv.sum()
print(f'mean NSTI {nsti.metadata_NSTI.mean():.3f}  median {nsti.metadata_NSTI.median():.3f}')
print(f'ASVs dropped at NSTI>{max_nsti}: {len(dropped)}/{len(nsti)}  reads dropped: {reads_dropped_frac:.1%}')
