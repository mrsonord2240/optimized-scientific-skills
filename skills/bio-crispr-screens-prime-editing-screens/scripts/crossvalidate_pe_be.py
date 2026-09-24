# Purpose: intersect PE and BE screen hits on the intended variant and flag high-confidence ones
#          (both FDR < cutoff and the same direction of effect).
# Inputs:  two TSVs, each with variant_id plus <prefix>_fdr and <prefix>_lfc columns
#          (be_fdr/be_lfc in the BE file, pe_fdr/pe_lfc in the PE file).
# Output:  merged TSV with a high_confidence column; prints the counts.
# Usage:   python crossvalidate_pe_be.py be_screen_hits.tsv pe_screen_hits.tsv [--fdr 0.05] [--out pe_be_concordance.tsv]
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('be_tsv', help='BE screen output (target conversion + bystander)')
ap.add_argument('pe_tsv', help='PE screen output (intended edit + scaffold-incorp + indel)')
ap.add_argument('--fdr', type=float, default=0.05)
ap.add_argument('--out', default='pe_be_concordance.tsv')
args = ap.parse_args()

be_hits = pd.read_csv(args.be_tsv, sep='\t')
pe_hits = pd.read_csv(args.pe_tsv, sep='\t')

# Intersect on intended variant
concordant = be_hits.merge(pe_hits, on='variant_id', suffixes=('_be', '_pe'))
# Filter to high-confidence: both methods call variant + same direction
concordant['high_confidence'] = ((concordant['be_fdr'] < args.fdr) & (concordant['pe_fdr'] < args.fdr) &
                                 (np.sign(concordant['be_lfc']) == np.sign(concordant['pe_lfc'])))
concordant.to_csv(args.out, sep='\t', index=False)
print(f'{len(concordant)} shared variants, {int(concordant["high_confidence"].sum())} high-confidence; wrote {args.out}')
